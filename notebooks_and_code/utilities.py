import phidl.geometry as pg
from phidl.quickplotter import quickplot2 as qp2
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from phidl.path import Path,CrossSection
import phidl.path as pth
import numpy as np
import shapely

from centerline.geometry import Centerline
from scipy.spatial import cKDTree
import warnings
import networkx as nx
import time
def phidl_to_shapely(dev):
    pols=dev.get_polygons()
    sh_pol = shapely.Polygon(pols[0])
    if len(pols)==1:
        return sh_pol
    else:
        for p in pols[1:]:
            sh_pol = sh_pol.union(shapely.Polygon(p))
        return sh_pol
def shapely_to_phidl(pol_round):
        rounded=pg.Device()
        rounded.add_polygon( np.asarray(pol_round.exterior.coords))
        if len(pol_round.interiors)>0:
            interior=pg.Device()
            for hole in pol_round.interiors:
                interior.add_polygon(hole.coords)
            rounded=pg.boolean(rounded,interior,operation='A-B')
        return rounded
def multilinestrings_to_phidl(multilinestring,layer=3,width=.1):
    dev=pg.Device()
    x=pth.CrossSection()
    x.add(width=width,layer=layer)
    for l in multilinestring.geoms:
        dev<<pth.Path(l.coords).extrude(x)
    return dev
def merge_line(line):
    new_linestring=[]
    for i in range(len(line)-1):
        l1=line[i]
        new_linestring.append(shapely.LineString(l1))
        l2=line[i+1]
        l1_a,l1_b,l2_a,l2_b=l1[0,:],l1[-1,:],l2[0,:],l2[-1,:]
        connect_table=np.sum(np.array([l1_a-l2_a,l1_a-l2_b,l1_b-l2_a,l1_b-l2_b])**2,axis=1) #want l1_b-l2_a to be min
        if np.min(connect_table)!=0: #there's a gap between lines
            if np.argmin(connect_table)==0:
                new_linestring.append(shapely.LineString([shapely.Point(l1_a),shapely.Point(l2_a)]))
            if np.argmin(connect_table)==1:
                new_linestring.append(shapely.LineString([shapely.Point(l1_a),shapely.Point(l2_b)]))
            if np.argmin(connect_table)==2:
                new_linestring.append(shapely.LineString([shapely.Point(l1_b),shapely.Point(l2_a)]))
            if np.argmin(connect_table)==3:
                new_linestring.append(shapely.LineString([shapely.Point(l1_b),shapely.Point(l2_b)]))
    new_linestring.append(shapely.LineString(line[-1]))
    return shapely.line_merge(shapely.unary_union(new_linestring))
def smooth_route(route,resample=5,num_final_coords=1000):
    from geomdl import NURBS,knotvector
    # Create a 3-dimensional B-spline Curve
    curve = NURBS.Curve()

    # Set degree
    curve.degree = 3
    curve.ctrlpts = np.asarray(shapely.LineString(route.interpolate(np.linspace(0,1,int(route.length/resample)),normalized=True)).coords)
    curve.knotvector = knotvector.generate(curve.degree,len(curve.ctrlpts))
    # Set knot vector

    # Set evaluation delta (controls the number of curve points)
    curve.delta = 1/num_final_coords

    # Get curve points (the curve will be automatically evaluated)
    curve_points = np.asarray(curve.evalpts)[:,:]
    return shapely.LineString(curve_points)

class RoutingLayout():#todo: make the loops with normals by using the distance function to check closest border to split lines equally max(linesgap, distance/n_lines)
    def __init__(self,device,border_tolerance,skeleton_fuse_dist=20,override_auto_skel_dist=False,skel_interp_dist=1,remove_dust_thresh=1e-1):
        assert type(device)==shapely.geometry.polygon.Polygon
        self.device=device
        self.shaved_device=device.buffer(-border_tolerance)

        #self.place_pads(pads) #extrude pads from shaved layout

        #create skeleton of the polygon (midlines)
        size=np.median(np.sqrt(np.sum(np.diff(np.asarray(self.shaved_device.envelope.boundary.coords),axis=0)**2,axis=1)))
        interp_dist = size/500
        if override_auto_skel_dist:
            interp_dist=skel_interp_dist
        cline=Centerline(self.shaved_device,interpolation_distance=interp_dist)
        cline=shapely.ops.linemerge(cline.geometry)
        self.skeleton=shapely.ops.linemerge(self.remove_dangling_lines(cline))
        self.skeleton=shapely.ops.linemerge(self.fuse_short_skeleton_segment(distance=skeleton_fuse_dist)) #sometimes two very close waypoints are created, which are topologically not necessary. This function snaps them and interpolates the incoming skeletons together
        self.skeleton=shapely.ops.linemerge([el for el in self.skeleton.geoms if el.length>remove_dust_thresh])
        self.electrodes_xy=[]
        self.pads_xy=[]
    def remove_dangling_lines(self,obj):
        '''self.obj is a shapely multiline'''
        xy_s=[]
        for line in obj.geoms:
            for g in line.boundary.geoms:
                p=np.asarray(g.coords).T
                xy_s.append(p)
        xy_s=np.asarray(xy_s)[:,:,0] #kill last dimension
        tree=cKDTree(xy_s)
        points_with_neighbors=np.unique(np.concatenate(list(tree.query_pairs(r=1e-5)))) 
        points_without_neighbors = [el//2 for el in range(2*len(obj.geoms)) if el not in points_with_neighbors]
        #we collect 2 points for each line edges in the kdtree, so line indices (who dangle) are obtained by //2
        new_skel=[obj.geoms[el] for el in range(len(obj.geoms)) if el not in points_without_neighbors]
        return shapely.ops.linemerge(shapely.MultiLineString(new_skel))
    
    #get waypoint to check closest border to understand what radius to look at. Should it do a square instead?
    
    def _lay_down_lines(self,d,n):
        collect=[]
        #d=self.line_to_line_distance
        #n=self.max_parallel_lines
        midline_polygon=self.skeleton.buffer(.1) #turn it into a polygon
        t1=time.time()
        if n%2==0: #odd line number
            gap_sizes = np.arange(d/2.,d/2.+n/2*d+1,d) #both sides of midline get only half a gap, since midline will be removed
        else:
            gap_sizes = np.arange(d,int(n/2.)*d+1,d)
        
        for i in gap_sizes: #start distance d from mid line and buffer down
            for pol in midline_polygon.interiors:
                route_line=shapely.Polygon(pol).buffer(-i).buffer(-5).buffer(5) #Make line to polygon, buffer it down, grab its outline
                if (route_line.geom_type!='MultiPolygon'):
                    if self.shaved_device.contains(route_line.exterior):
                        collect.append(np.asarray(route_line.exterior.coords))
                    else:
                        left_over_line=route_line.exterior.intersection(self.shaved_device)
                        xys=self.get_coords(left_over_line)
                        if type(xys)==list:
                            for coords in xys:
                                collect.append(coords)
                        elif len(xys)>0:
                            collect.append(xys)

                #     collect.append(np.asarray( shapelypolygon.difference(shapely.Polygon(route_line)).exterior.coords))
            route_line=shapely.Polygon(midline_polygon.exterior).buffer(i).buffer(5).buffer(-5).exterior #Make line to polygon, buffer it down, grab its outline
            if self.shaved_device.contains(route_line):
                collect.append(np.asarray(route_line.coords))
            else:
                left_over_line=route_line.intersection(self.shaved_device)
                xys=self.get_coords(left_over_line)
                if type(xys)==list:
                    for coords in xys:
                        collect.append(coords)
                elif len(xys)>0:
                    collect.append(xys)
        
        t2=time.time()
        self.possible_routes=shapely.MultiLineString([shapely.LineString(el) for el in collect])
        if n%2!=0:
            self.possible_routes = self.skeleton.union(self.possible_routes) #put the midline as a route.
        self.make_grids_at_waypoints() #at every crossing, rotate the local lines to make an exchanger grid.
        t3=time.time()
        self.lines_to_add_for_crisscross=shapely.ops.linemerge(shapely.ops.unary_union(self.lines_to_add_for_crisscross)) #simplify topo
        diff_ab=self.possible_routes.difference(self.lines_to_add_for_crisscross) #chop routes by crosscross to create separate line bits
        diff_ba=self.lines_to_add_for_crisscross.difference(self.possible_routes) #chop criss cross by routes to create separate line bits
        self.all_route_lines=shapely.ops.linemerge(diff_ab.union(diff_ba)) #collect all lines in a multistring
        t4=time.time()
        print(t2-t1,t3-t2,t4-t3)
    def lay_down_lines(self,n_lines,gap_um,adaptive=False,add_border=True):
        def arclength_line_and_normal(linestring,resolution=5):
            '''resolution in microns, float'''
            arclength = np.linspace(0,1, int(linestring.length/resolution))
            arcline=shapely.LineString(linestring.interpolate(arclength,normalized=True))
            xy=np.asarray(arcline.xy).T
            norm_xy = np.gradient(xy,axis=0)
            norm_xy = norm_xy/np.sqrt(np.sum(norm_xy**2,axis=1)).reshape(-1,1)
            return xy,np.vstack((-norm_xy[:,1],norm_xy[:,0])).T
        def clean_up_nonsimple(xy,norm,distance,buffer_radius=5):
            '''removes hairpins and buffers out dirty corners'''
            pts=xy
            pts_in = xy-norm*distance_threshold.reshape(-1,1)
            #internal loop
            #linemerge puts together contiguous non intersection line bits. these are long. the artefacts create hairpins, which are short
            #we remove shortest line elements until the linemerge creates a singe linestring => no self-interesections    

            all_linestrings=shapely.line_merge(shapely.unary_union([shapely.LineString([p1,p2]) for p1,p2 in zip(pts_in,np.roll(pts_in,axis=0,shift=1))]))
            if all_linestrings.geom_type=='LineString':
                pruned_ring = shapely.LinearRing(all_linestrings)
                pruned_ring=pruned_ring.buffer(buffer_radius+1e-3).buffer(-buffer_radius).interiors[0]
                if not pruned_ring.is_ccw:
                    pruned_ring=pruned_ring.reverse()
            else:
                lengths=[el.length for el in all_linestrings.geoms]
                ordered_lengths=np.argsort(lengths)
                ctr=1
                pruned_linestrings=all_linestrings
                while pruned_linestrings.geom_type=='MultiLineString':
                    pruned_linestrings=[all_linestrings.geoms[el] for el in range(len(lengths)) if el in ordered_lengths[ctr:]]  #remove shortest artefacts one by one
                    pruned_linestrings=shapely.line_merge(shapely.unary_union(pruned_linestrings))
                    ctr+=1
                pruned_ring = shapely.LinearRing(pruned_linestrings)
                pruned_ring=pruned_ring.buffer(buffer_radius+1e-3).buffer(-buffer_radius).interiors[0]
                if not pruned_ring.is_ccw:
                    pruned_ring=pruned_ring.reverse()
    
            return pruned_ring
        self.simp_skel = self.skeleton.simplify(.1).segmentize(1) #speeds up distance computation.
        self.simp_skel_buff = self.simp_skel.buffer(.1) #make it a polygon
        boundary_growths = n_lines//2
        self.possible_routes=[]
        all_lines=[]
        if adaptive:
            for loopidx,intloop in enumerate(self.simp_skel_buff.boundary.geoms):
                linestring=shapely.Polygon(intloop).buffer(2*gap_um).buffer(-4*gap_um).buffer(2*gap_um).exterior #round corners to make sure normal curve is smooth
                if not linestring.is_ccw:
                    linestring=linestring.reverse()
                sign=-1
                if loopidx==0:
                    sign=1 #orient the normals outward when taking the exterior of skelbuff.
                xy,norm=arclength_line_and_normal(linestring,resolution=2)
                norm*=sign
                shaved_outline=shapely.unary_union(self.shaved_device.boundary)
                lines=[]
                for i in range(boundary_growths-1):
                    distance = np.asarray([shapely.Point(x,y).distance(shaved_outline) for x,y in xy])
                    distance_proposed=distance/(boundary_growths-i)
                    distance_threshold=np.copy(distance_proposed)
                    distance_threshold[distance_threshold<gap_um]=gap_um
            #distance_threshold[distance_threshold<(2*gap_um)]=1.1*distance[distance_threshold<(2*gap_um)] #if it cant fit all curves, divide equally for fewer curves
                    factor = distance/distance_threshold #smooth out this factor such that loop curve more closely follows device boundary
            #smoothed_factor=medfilt(factor,31)
            #smoothed_distance = distance/smoothed_factor
                    newline = clean_up_nonsimple(xy,norm,distance_threshold,buffer_radius=5)
                    lines.append(newline)
           # plt.plot(distance)
                    xy,norm=arclength_line_and_normal(newline,resolution=2)
                    norm*=sign
        
                device_loop = shapely.Polygon([el for el in self.shaved_device.interiors if shapely.Polygon(linestring).contains(el)][0])
                lines=shapely.unary_union(lines)
            #lines=lines.difference(device_loop.buffer(gap_um)) #kill lines too close to device loop, whose exterior is also a electrode line
            #lines=lines.intersection(sh.Polygon(linestring).buffer(-gap_um)) #kill lines too close to skeleton midline
                all_lines.append(lines)
    
        if not adaptive:
            for i in range(boundary_growths):
                lines=self.simp_skel_buff.buffer((i+1)*gap_um).boundary.intersection(self.shaved_device)
                all_lines.append(shapely.line_merge(lines))
                
            
        all_lines.append(self.simp_skel)
        self.possible_routes=shapely.unary_union(all_lines) #dont prune to make longer waypoint cross-curves
        self.make_grids_at_waypoints(waypoint_size=1)
        all_lines=shapely.unary_union(all_lines).intersection(self.shaved_device.buffer(-gap_um*.9)) #clean up edges              
        if add_border:
        	all_lines=shapely.unary_union(all_lines.union(self.shaved_device.boundary))
        all_lines=shapely.unary_union(all_lines.union(shapely.unary_union(self.lines_to_add_for_crisscross)))
        self.all_route_lines=shapely.line_merge(all_lines)
    def add_electrodes(self,list_of_boundaries,list_of_points_to_connect_center=[],clearance=5,type='electrode'):
        #clearance is in microns
        if type not in ['electrode','pad']:
            warnings.warn('Incorrect element: should one of [electrode, pad]')
        cutouts=[]
        all_new_paths=[]
        centers=[]
        for boundary,point_list in zip(list_of_boundaries,list_of_points_to_connect_center):
            object_clearance = shapely.Polygon(boundary).buffer(clearance)
            cutouts.append(object_clearance)
            outline_clearance=boundary.buffer(clearance).exterior
            center_xy = list(boundary.centroid.coords)[0]
            centers.append(center_xy)
            attachment_to_routes = self.all_route_lines.intersection(object_clearance) #the route lines covered by the object
    #this attaches all intersections between object clearance boundary to its center line, to route paths through electrodes/pads.
            new_paths=shapely.MultiLineString([shapely.LineString( [ shapely.Point(center_xy), el]) for el in attachment_to_routes.boundary.geoms])
            if point_list.geom_type=='MultiPoint':
                new_paths=new_paths.union([shapely.LineString( [shapely.Point(center_xy),el]) for el in point_list.geoms])
            all_new_paths.append(new_paths)
        self.all_route_lines=self.all_route_lines.difference(shapely.unary_union(cutouts))
        self.all_route_lines=shapely.unary_union(self.all_route_lines.union(shapely.unary_union(all_new_paths)))
        if type=='electrode':
            self.electrodes_xy=centers
        if type=='pad':
            self.pads_xy=centers #important for graph path
    def visualize_router_qp(self,width=1,layer=3):
        dev=pg.Device()
        x=pth.CrossSection()
        x.add(width=width,layer=layer)
        for l in self.all_route_lines.geoms:
            dev<<pth.Path(l.coords).extrude(x)
        return dev
    def route(self,resample_for_smooth=10):
        line_end_coords=[]
        self.all_lines=[]
        for l in self.all_route_lines.geoms:
                bdy_coords=[list(el.coords)[0] for el in l.boundary.geoms]
                if len(bdy_coords)>0:
                    line_end_coords.append(bdy_coords) #end point coords of lines
                    self.all_lines.append(l) #collect all lines into one list
        #there are faulty lines in the merge, which fold back onto themselves (first point = last point), which technically does not have a boundary and is a closed polygon
        #eliminate those
        self.end_pts=np.concatenate(line_end_coords)
        end_to_end_connections= np.asarray([(2*i,2*i+1) for i in range(len(self.end_pts)//2)]) #every consecutive 2 points are connected together by the LineString
        
        tree1=cKDTree(self.end_pts) #the end and beginning point of two consecutive lines is considered as 2 distinct points. We do a renaming here for the purpose of the graph so these appear as the same node.
        nodes_to_merge=tree1.query_ball_tree(tree1,r=1e-3)
        for nodes in nodes_to_merge:
            if len(nodes)>0:
                for replaced in nodes[1:]: #replace all other neighbor nodes by the first node index, so intersections dont count as multiple
                    end_to_end_connections[end_to_end_connections==replaced]=nodes[0]
        
        self.end_to_end_connections=end_to_end_connections #renamed
        self_connections=np.where(self.end_to_end_connections[:,0]-self.end_to_end_connections[:,1] == 0)[0]
        self.end_to_end_connections = np.asarray([self.end_to_end_connections[i] for i in range(len(end_to_end_connections)) if i not in self_connections])
        
        self.all_lines=shapely.MultiLineString([self.all_lines[i] for i in range(len(self.all_lines)) if i not in self_connections])
        self.G=nx.Graph()
        self.G.add_edges_from(self.end_to_end_connections)
        danglers = [k for k,v in self.G.degree if v==1]
        while len(danglers)>0:
            self.G.remove_nodes_from([k for k,v in self.G.degree if v==1])
            danglers = [k for k,v in self.G.degree if v==1]
            
        source,drain=int(1e15),int(1e16)
        all_elec_nodes=cKDTree(np.asarray(self.electrodes_xy)).query_ball_tree(cKDTree(self.end_pts),r=1e-3) #the center node is degenerate, unique index is collapsed to a single one inside graph
        all_elec_nodes_clean=[[el for el in degenerate_elec_index if el in self.G.nodes][0] for degenerate_elec_index in all_elec_nodes]
        all_elec_edges = [(el,source) for el in all_elec_nodes_clean]
    
        all_pad_nodes=cKDTree(np.asarray(self.pads_xy)).query_ball_tree(cKDTree(self.end_pts),r=1e-3) 
        all_pad_nodes_clean=[[el for el in degenerate_pad_index if el in self.G.nodes][0] for degenerate_pad_index in all_pad_nodes]
        all_pad_edges = [(el,drain) for el in all_pad_nodes_clean]
        self.G.add_edges_from(all_pad_edges)
        self.G.add_edges_from(all_elec_edges)
        
        
        paths=[el[1:-1] for el in list(nx.node_disjoint_paths(self.G,1e15,1e16))]
        if len(paths)==len(self.electrodes_xy):
            print("Routing successful")
        self.all_elec_pad_routes=[]
        self.all_elec_pad_routes_smooth=[]
        #shortest_paths=[]
        
        #for idx,p in enumerate(paths):
        #    other_paths=[p for pidx,p in enumerate(paths) if pidx!=idx]
        #    subG=self.G.copy()
        #    for op in other_paths:
        #        subG.remove_nodes_from(op)
        #    shortest_paths.append(nx.shortest_path(subG,p[0],p[-1]))
        
        for p in paths:
            linebits=[]
            for i in range(len(p)-1):
                line_idx=np.where(np.sum(self.end_to_end_connections==p[i],axis=1) + np.sum(self.end_to_end_connections==p[i+1],axis=1) == 2)[0][0]
                linebits.append(self.all_lines.geoms[line_idx])
            linebits=merge_line([np.asarray(l.coords) for l in linebits])
            self.all_elec_pad_routes_smooth.append(smooth_route(linebits,resample=resample_for_smooth,num_final_coords=2*int(linebits.length)))
            self.all_elec_pad_routes.append(linebits)
        self.paths=paths
        #self.shortest_paths=shortest_paths
    def plot_possible_routes(self):
        for el in self.possible_routes.geoms:
            x,y=np.asarray(el.coords).T
            plt.plot(x,y,color='k')
    def get_coords(self,obj):
        if hasattr(obj,'geoms'):
            coords=[]
            for g in obj.geoms:
                coords.append(np.asarray(g.coords))
        else:
            coords=np.asarray(obj.coords)
        return coords
    def make_grids_at_waypoints(self,waypoint_size=1,rotation_angle=180):
        #All biffurcations in the skeleton indicate points at which routes need to jump to different lines.
        #to do that, we need to provide a local grid. This is done by cutting out all the lines around
        #a biffurcation point, and rotating them in place and adding them to the design. This creates
        #a criss-cross pattern.
        self.borders=shapely.MultiLineString([self.device.exterior]+[el for el in self.device.interiors])
        #waypoints=self.skeleton.boundary #gets all end points of line segments. Check if not double counted
        #waypoints are sometimes double counted or too close and are not considered as a single one. Clean it up of closer than 1e-3
        waypoints = shapely.unary_union([g.boundary for g in self.skeleton.geoms])
        tree=cKDTree((np.asarray([list(p.coords)[0] for p in waypoints.geoms])))
        exclude=[v for k,v in tree.query_pairs(1e-3)] #pick second waypoint and exlcude it (closer than tolerance to first)
        self.waypoints=shapely.MultiPoint([g for idx,g in enumerate(waypoints.geoms) if idx not in exclude])    
    
        self.lines_to_add_for_crisscross = []
        for wp in self.waypoints.geoms:
            number_of_branches=np.sum([el.intersects(wp.buffer(waypoint_size)) for el in self.skeleton.geoms])
            rotation_angle=180/number_of_branches
            x,y=np.asarray(wp.coords).T
            radius_to_border = wp.distance(self.borders)*1.05
            lines_around_wp = wp.buffer(radius_to_border).intersection(self.possible_routes) #make a circle around waypoint and collect all lines in it
            rotated=shapely.affinity.rotate(lines_around_wp,rotation_angle,origin=(x,y))
            #rotated_clean=self.remove_dangling_lines(rotated.union(lines_around_wp)) #remove dangler lines in the cross cross
            #rotated_clean=rotated.intersection(rotated_clean) #retrieve original extra lines, to not carry a copy of original 'lines around wp' lines inside the route lines
            #rotated_clean=unary_union([el for el in rotated_clean.geoms if el.geom_type=='LineString']) #points are not necessary here 
            #self.lines_to_add_for_crisscross.append(rotated_clean)#reduce number of segments if contiguous
            self.lines_to_add_for_crisscross.append(rotated)
        self.lines_to_add_for_crisscross = self.shaved_device.intersection(self.lines_to_add_for_crisscross)

        #clean up these lines. boolean cut them, remove dangles, etc.
    def fuse_short_skeleton_segment(self,distance=20):
        '''Sometimes the skeleton has short inner segments that can be collapsed. This function does that.
        It is important to make cleaner waypoints for routing, as a short segment inside the skeleton produces
        two waypoints that are closed to each other >--<'''
        def stretch_to_point(curve,point):
            closest_curve_endp=np.argmin([point.distance(el) for el in curve.boundary.geoms])
            stretch_end = curve.boundary.geoms[closest_curve_endp]
            anchor = curve.boundary.geoms[int( np.abs(1-closest_curve_endp))]
            
            curve_vector = np.asarray(stretch_end.coords)-np.asarray(anchor.coords) 
            stretched_vector = np.asarray(point.coords)-np.asarray(anchor.coords) 
            
            curve_r = np.sqrt(np.sum(curve_vector**2))
            stretched_r = np.sqrt(np.sum(stretched_vector**2))
            
            curve_angle = np.arctan2(*curve_vector[0,:][::-1]) #invert order sequence (arctan(y,x))
            stretched_angle = np.arctan2(*stretched_vector[0,:][::-1]) #invert order sequence (arctan(y,x))
            
            angle_diff = stretched_angle-curve_angle
            
            stretched_curve=shapely.affinity.scale(curve,xfact=stretched_r/curve_r, yfact=stretched_r/curve_r, origin=tuple(np.asarray(anchor.coords)[0]))
            rotated_stretched_curve=shapely.affinity.rotate(stretched_curve,angle=angle_diff*180/np.pi, origin = tuple(np.asarray(anchor.coords)[0]))
        
            #stretching it like that moves the curve off skeleton midline, and the routing loops may collide with device borders.
            a=curve
            b=rotated_stretched_curve
        
            pa,pb=np.asarray(a.coords),np.asarray(b.coords)
            if closest_curve_endp==0:
                pa,pb=pa[::-1,:],pb[::-1,:]
            t=np.linspace(0,1,len(pa))

            switch_curve_at_fraction = np.max( (0,1-(distance//5)/b.length)) #start fusing to point distance microns before target point
            fac=0.5*(np.tanh( (t-switch_curve_at_fraction)*30)+1).reshape(len(t),1)
            
            fac-=fac[0]
            fac/=fac[-1]
            
            c = pa*(1-fac)+pb*fac
            
            return shapely.LineString(c)
        boundaries=[]
        skeleton_lines=[el for el in self.skeleton.geoms]
        for g in self.skeleton.geoms:
            for b in g.boundary.geoms:
                boundaries.append(b.xy)
        edges=cKDTree(np.asarray(boundaries)[:,:,0]).query_pairs(distance)

        G=nx.Graph()
        G.add_edges_from(edges)
        clusters=list(nx.connected_components(G))
        throw_idx=[]#remove internal danglers
        for c in clusters:
            mean_px,mean_py = np.mean(np.asarray([boundaries[e] for e in c])[:,:,0],axis=0)
            line_idx=np.unique([e//2 for e in c])
            for lidx in line_idx:
                line=skeleton_lines[lidx]
                if line.length>distance:
                    arranged=stretch_to_point(line,shapely.Point(mean_px,mean_py))
                    skeleton_lines[lidx]=arranged
                else:
                    throw_idx.append(lidx)
        new_skel = [el for idx,el in enumerate(skeleton_lines) if idx not in throw_idx]
        return shapely.line_merge(shapely.unary_union(new_skel))

    def place_pads(self,phidl_geometry):
        pad_pol = []
        for p in phidl_geometry.get_polygons():
            pad_pol.append(shapely.Polygon(p))
        self.shaved_device=self.shaved_device.difference(shapely.ops.unary_union(pad_pol))
