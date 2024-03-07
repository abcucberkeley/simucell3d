#include "ball_pivoting_algorithm.hpp"


/*
    Given a poisson point cloud, this class reconstructs a triangulated surface by using the 
    ball pivoting algorithm. The theory behind this algorithm can be found in the original paper:
    "The ball-pivoting algorithm for surface reconstruction", by Bernardini et al, DOI: 10.1109/2945.817351
*/


//----------------------------------------------------------------------------------------------
ball_pivoting_algorithm::ball_pivoting_algorithm(
        const std::vector<oriented_point>& poisson_point_cloud, 
        const double l_min, 
        const double min_x, const double min_y, const double min_z, 
        const double max_x, const double max_y, const double max_z
    ) noexcept(false): ball_radius_(1.7 * l_min), node_lst_(poisson_point_cloud){
    ball_radius_squared_ = ball_radius_*ball_radius_;

    //Run some checks
    assert(ball_radius_ > 0.);
    assert(node_lst_.size() >= 4);
    assert(node_lst_.size()  < std::numeric_limits<unsigned>::max());

    //Create the space partitionning grid
    grid_ = uspg_4d<oriented_point>(
        min_x - ball_radius_ * 2., // Without this padding, the grid will not be able to find the points at the border of the grid
        min_y - ball_radius_ * 2., //near a circumspherecenter outside the grid.
        min_z - ball_radius_ * 2., 
        max_x + ball_radius_ * 2., 
        max_y + ball_radius_ * 2., 
        max_z + ball_radius_ * 2., 
        ball_radius_, node_lst_.size());

    //Populate the grid and set the point ids
    for(unsigned point_id = 0; point_id < node_lst_.size(); point_id ++){

        //Get a ref to the point object
        oriented_point& p = node_lst_[point_id];

        //Set the id of the point
        p.id_ = point_id;

        //And place the point in the grid
        grid_.place_object(p, p.position_);
    }
    
}
//----------------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------------
//Find a triangle from which the surface expansion can be started. 
//Return false if it cannot find a seed triangle
void ball_pivoting_algorithm::find_seed_triangle() noexcept(false){

    //Basically this methods loop over the point cloud, and try to connect each point with
    //2 other neighboring points to form a triangle. As soon as it is done, it returns the first created face

    //Create a vector of node indices
    std::vector<unsigned> index_lst(node_lst_.size());
    std::iota(index_lst.begin(), index_lst.end(), 0);

    //Shuffle the indices, otherwise the seed triangle will be always the same for a given point cloud
    std::random_shuffle(index_lst.begin(), index_lst.end());

    //Loop over the points of the poisson point cloud
    for(const auto node_index: index_lst){
        auto& a = node_lst_[node_index];

        //Get the points nearby
        const auto nearby_point_lst = grid_.get_neighborhood(a.position_);

        for(auto b: nearby_point_lst){
        for(auto c: nearby_point_lst){

            //Make sure that all the points are distinct in the triplet of points
            if(a.id_ == b.id_ || a.id_ == c.id_ || b.id_ == c.id_) continue;

            //Make sure this points have normals pointing in the same direction
            if(!ball_pivoting_algorithm::points_coherently_oriented(a, b, c)) continue;

            //Get the center of the sphere of radius = ball_radius_ passing by the 3 points
            //If this method doesn't return anything than this circumsphere does not exist
            auto circum_sphere_center_opt = compute_circum_sphere_center(a, b, c);

            //If a circumsphere with the correct radius could not be formed
            if(!circum_sphere_center_opt) continue;
            const vec3 circum_sphere_center = circum_sphere_center_opt.value();

            //Make sure there are no other nodes inside of the circumsphere
            if(!circum_sphere_is_empty(circum_sphere_center, a, b, c)) continue;

            //If all the checks have been passed with success, create 3 new edges
            auto e1_id = get_edge(a, b);   assert(e1_id < edge_lst_.size());
            auto e2_id = get_edge(b, c);   assert(e2_id < edge_lst_.size());
            auto e3_id = get_edge(c, a);   assert(e3_id < edge_lst_.size());

            edge& e1 = edge_lst_[e1_id];
            edge& e2 = edge_lst_[e2_id];
            edge& e3 = edge_lst_[e3_id];

            

            //Create a face with the 3 nodes and make sure its normal is correctly oriented
            const vec3 correct_normal = ball_pivoting_algorithm::compute_oriented_normal(a, b, c);
            const vec3 initial_normal = (a.position_ - b.position_).cross(a.position_ - c.position_);

            //Create the face with the correct orrientation            
            face& f = (initial_normal.dot(correct_normal) < 0) ? face_lst_.emplace_back(a.id_, c.id_, b.id_, 0) : face_lst_.emplace_back(a.id_, b.id_, c.id_, 0); 

            //Connect the edges with the newly created face
            e1.add_face(f.get_local_id());
            e2.add_face(f.get_local_id());
            e3.add_face(f.get_local_id());

            //Connect the nodes with this face
            node_face_lst_[a.id_ ].push_back(f.get_local_id());
            node_face_lst_[b.id_ ].push_back(f.get_local_id());
            node_face_lst_[c.id_ ].push_back(f.get_local_id());

            //Indicate to the nodes that they are being used by these edges
            node_edge_lst_[a.id_].push_back(e1_id);
            node_edge_lst_[a.id_].push_back(e3_id);

            node_edge_lst_[b.id_].push_back(e1_id);
            node_edge_lst_[b.id_].push_back(e2_id);

            node_edge_lst_[c.id_].push_back(e2_id);
            node_edge_lst_[c.id_].push_back(e3_id);

            //Store the 3 newly created edges into the edge front
            edge_id_front_lst_.push_back(e1_id);
            edge_id_front_lst_.push_back(e2_id);
            edge_id_front_lst_.push_back(e3_id);

            //Return that a seed triangle has been found
            return;
        }}        
    }

    //Throw an exception if no seed triangle could be found
    throw bpa_exception("No seed triangle found");
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Expand the triangulation from the seed triangle
void ball_pivoting_algorithm::expand_triangulation() noexcept{

    //If you raised this assertion error it means that you called this method without 
    //first running the find_seed_triangle() method
    assert(!edge_id_front_lst_.empty());

    while(!edge_id_front_lst_.empty()){

        //std::cout << "Edge front size: " << edge_id_front_lst_.size() << std::endl;

        //Pop the last edge of the edge_front
        auto edge_id = edge_id_front_lst_.back(); assert(edge_id < edge_lst_.size());
        edge_id_front_lst_.pop_back();
        edge& e = edge_lst_[edge_id];

        //If the edge has already been transformed into an inner edge
        if(e.is_manifold()) continue;

        //Find the best candidate node to form a new triangle with this edge
        std::optional<oriented_point> point_c_opt = find_candidate_node(e); 
        
        //If no good candidate was found to create a triangle with the edge e
        if(!point_c_opt) continue; //Go to the next edge
        oriented_point& point_c = point_c_opt.value();

        //Create the new face made of the nodes of the edge and the candidate point
        auto [point_a_id, point_b_id] = e.get_node_ids();
        assert(point_a_id < node_lst_.size() && point_b_id < node_lst_.size());
        oriented_point& point_a = node_lst_[point_a_id];
        oriented_point& point_b = node_lst_[point_b_id];

        //Update the edges of face_b
        auto e1_id = get_edge(point_a, point_b); assert(e1_id < edge_lst_.size()); 
        auto e2_id = get_edge(point_b, point_c); assert(e2_id < edge_lst_.size());
        auto e3_id = get_edge(point_c, point_a); assert(e3_id < edge_lst_.size()); 

        
        edge& e1 = edge_lst_[e1_id]; edge& e2 = edge_lst_[e2_id]; edge& e3 = edge_lst_[e3_id];

        //WARNING: This line should be checked if there are problems with the mesh
        if(e1.is_manifold() || e2.is_manifold() || e3.is_manifold()) continue;


        //Create a face with the 3 nodes and make sure its normal is correctly oriented
        const vec3 correct_normal = ball_pivoting_algorithm::compute_oriented_normal(point_c, point_a, point_b);
        const vec3 initial_normal = (point_a.position_ - point_c.position_).cross(point_b.position_ - point_c.position_);

        //Get the face_id
        const unsigned face_id = static_cast<unsigned>(face_lst_.size());

        //Create the face with the right node order such that its normal points outward            
        auto face_b = (initial_normal.dot(correct_normal) < 0.) ? 
        face_lst_.emplace_back(point_c.id_, point_b.id_, point_a.id_, face_id) : face_lst_.emplace_back(point_c.id_, point_a.id_, point_b.id_, face_id);
            


        //Connect the edges with the newly created face
        e1.add_face(face_id);
        e2.add_face(face_id);
        e3.add_face(face_id);

        //Connect the nodes with this face
        node_face_lst_[point_a.id_].push_back(face_b.get_local_id());
        node_face_lst_[point_b.id_].push_back(face_b.get_local_id());
        node_face_lst_[point_c.id_].push_back(face_b.get_local_id());


        node_edge_lst_[point_a.id_].push_back(e1_id);
        node_edge_lst_[point_a.id_].push_back(e3_id);

        node_edge_lst_[point_b.id_].push_back(e1_id);
        node_edge_lst_[point_b.id_].push_back(e2_id);

        node_edge_lst_[point_c.id_].push_back(e2_id);
        node_edge_lst_[point_c.id_].push_back(e3_id);


        //Store the 3 newly created edges into the edge front
        if(!e1.is_manifold()) edge_id_front_lst_.push_back(e1_id);
        if(!e2.is_manifold()) edge_id_front_lst_.push_back(e2_id);
        if(!e3.is_manifold()) edge_id_front_lst_.push_back(e3_id);
            
    }
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Check if the points have normals pointing in the same direction
bool ball_pivoting_algorithm::points_coherently_oriented(
    const oriented_point& a, 
    const oriented_point& b,
    const oriented_point& c
) noexcept{

    //Get the normal formed by the 3 nodes
    const vec3 normal =  ball_pivoting_algorithm::compute_oriented_normal(a,b,c);
    if(a.normal_.dot(normal) < 0.) return false;
    if(b.normal_.dot(normal) < 0.) return false;
    if(c.normal_.dot(normal) < 0.) return false;

    return true;
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Get the center of the sphere of radius = ball_radius_ passing by the 3 points
//If this method doesn't return anything then this circumsphere does not exist
std::optional<vec3> ball_pivoting_algorithm::compute_circum_sphere_center(
    const oriented_point& point_a, 
    const oriented_point& point_b,
    const oriented_point& point_c
) const noexcept{

    //Compute the length of the 3 edges of the triangle
    double a = (point_b.position_ - point_c.position_).norm();
    double b = (point_a.position_ - point_c.position_).norm();
    double c = (point_a.position_ - point_b.position_).norm();

    //Caculate the radius of the circumshpere
    const double squared_radius =   (std::pow(a,2) * std::pow(b,2) * std::pow(c,2)) / 
                                    ((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));


    //If the radius of the circumsphere exceeds the maximum ball radius
    if(squared_radius >= ball_radius_squared_) return std::nullopt;


    //Calculate the center of the circle positioned on the plane formed by the 3 nodes
    //and that at the same time overlaps with the 3 nodes.
    const vec3 ab = point_b.position_ - point_a.position_;
    const vec3 ac = point_c.position_ - point_a.position_;
    const vec3 normal = ab.cross(ac);

    const vec3 circle_circum_center_point =   point_a.position_ +
    ((normal.cross(ab) * ac.squared_norm()) + (ac.cross(normal) * ab.squared_norm())) / (2. * normal.squared_norm());


    //Now use the normal of the triangle to translate this point orthogonaly to the plane
    //formed by the triangle

    //Calculate by how much the circle circum center must be translated along the triangle normal
    const double translation_coeff = std::sqrt(ball_radius_squared_ - squared_radius);
    assert(std::isfinite(translation_coeff));


    //Get the normal of the triangle formed by this three nodes
    const vec3 oriented_normal = ball_pivoting_algorithm::compute_oriented_normal(point_a, point_b, point_c);


    //Compute the center of the circum shpere
    const vec3 sphere_circum_center_point = circle_circum_center_point + (oriented_normal * translation_coeff);
    
    
    return std::optional<vec3>{sphere_circum_center_point};
}


std::optional<vec3> ball_pivoting_algorithm::compute_circum_sphere_center(
    const face& f
) const noexcept{

    //Get the ids of the nodes making up the face
    auto [n1_id, n2_id, n3_id] = f.get_node_ids();

    assert(n1_id < node_lst_.size());
    assert(n2_id < node_lst_.size());
    assert(n3_id < node_lst_.size());

    const oriented_point& node_1 = node_lst_[n1_id];
    const oriented_point& node_2 = node_lst_[n2_id];
    const oriented_point& node_3 = node_lst_[n3_id];

    return compute_circum_sphere_center(node_1, node_2, node_3);
}

//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
//Compute the normal of the plane formed by the 3 points and make 
//that this normal is oriented in the same direction as the points' normals.
//The retruned normal is a unit vector
vec3 ball_pivoting_algorithm::compute_oriented_normal(
    const oriented_point& point_a, 
    const oriented_point& point_b,
    const oriented_point& point_c
) noexcept{

    //Compute the triangle of the face
    vec3 normal = (point_b.position_ - point_a.position_).cross(point_c.position_ - point_a.position_);

    //Compute the average normal of the nodes
    vec3 avg_normal = point_a.normal_.normalize() + point_b.normal_.normalize() + point_c.normal_.normalize();

    //Make sure the face normal has the same orientation as the average nodes normal
    if(normal.dot(avg_normal) < 0.) normal = normal * (-1.);

    return normal.normalize();
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Check if the circum_sphere with center circum_sphere_center and radius ball_radius_
//does not overlap with points that are point_a, point_b, or point_c..
bool ball_pivoting_algorithm::circum_sphere_is_empty(
        const vec3& circum_sphere_center, 
        const oriented_point& point_a, 
        const oriented_point& point_b, 
        const oriented_point& point_c
    ) const noexcept{

    //Get all the points that are close to the circum_sphere_center point
    auto nearby_point_lst = grid_.get_neighborhood(circum_sphere_center);

    //Loop over the nearby points
    for(oriented_point p: nearby_point_lst){

        //Check that they are not part of the given triplet
        if(p.id_ == point_a.id_ || p.id_ == point_b.id_ || p.id_ == point_c.id_ ) continue;

        //Compute the distance with the center of the circum sphere
        const double squared_distance = (p.position_ - circum_sphere_center).squared_norm();

        //If the point is in the circum sphere
        if(squared_distance <= ball_radius_squared_) return false;
    }

    return true;
}
//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
//Returns the id of the edge formed by 2 points.
unsigned ball_pivoting_algorithm::get_edge(const oriented_point& point_a, const oriented_point& point_b) noexcept{
    return get_edge(point_a.id_, point_b.id_);
}


unsigned ball_pivoting_algorithm::get_edge(const unsigned point_1_id, const unsigned point_2_id) noexcept{
    
    //Find if the edge has already been created
    auto edge_it = std::find(edge_lst_.begin(), edge_lst_.end(), edge(point_1_id, point_2_id));

    //If the edge doesn't exist yet we create it and add it to the edge_lst_
    if(edge_it == edge_lst_.end()){
        edge_lst_.emplace_back(point_1_id, point_2_id);
        return static_cast<unsigned>(edge_lst_.size()) -1;
    }

    return static_cast<unsigned>(std::distance(edge_lst_.begin(), edge_it));

}
//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
//Returns the id of the best node with which to connect the edge to form a triangle
std::optional<oriented_point> ball_pivoting_algorithm::find_candidate_node(const edge& edge) const noexcept{

    //Get the first face of the edge
    const auto face_a_id = edge.f1();
    assert(face_a_id < face_lst_.size());
    const face& face_a = face_lst_[face_a_id];

    //Get references to the nodes making up the edge
    auto [point_a_id, point_b_id] = edge.get_node_ids();
    assert(point_a_id < node_lst_.size() && point_b_id < node_lst_.size());
    
    const oriented_point& point_a = node_lst_[point_a_id];
    const oriented_point& point_b = node_lst_[point_b_id];

    //Get a vector of the edge
    const vec3 edge_vec = (point_b.position_ - point_a.position_);

    //Calculate the center of the edge
    const vec3 edge_center_vec = (edge_vec * 0.5) + point_a.position_;

    //Get the node opposite to the edge in the face_a
    auto point_c_id = face_a.get_opposite_node(point_a_id, point_b_id);
    assert(point_c_id < node_lst_.size());
    const oriented_point& point_c = node_lst_[point_c_id];

    //Compute the center of circum sphere of the face a 
    auto face_a_circum_center_op = compute_circum_sphere_center(face_a);
    assert(face_a_circum_center_op);
    const vec3 face_a_circum_center = face_a_circum_center_op.value();

    //The node id of best candidate node to form a triangle with the edge
    std::optional<oriented_point> best_candidate;

    //WARNING: This is set to M_PI in the first version
    double min_theta = std::numeric_limits<double>::infinity();

    //Loop over all the nodes located nearby the edge center
    const auto neighboring_nodes = grid_.get_neighborhood(edge_center_vec);
    for(const oriented_point p: neighboring_nodes){

        //Make sure the candidate point is not part of face a
        if(p.id_ == point_a_id || p.id_ == point_b_id || p.id_ == point_c_id) continue;

        //Make sure the node is not already inside the triangulated surface
        if(node_is_inner_vertex(p)) continue;

        //Check that the normals of the faces {point_a, point_b, point_c}, and 
        //{point_a, point_b, p} point in the same direction
        if(!ball_pivoting_algorithm::node_is_compatible_with_edge(p, point_a, point_b, point_c)) continue;

        //Compute the circumsphere center of the face {point_a, point_b, p}
        auto circum_sphere_center_opt = compute_circum_sphere_center(point_a, point_b, p);

        //If a circumsphere with a radius == ball_radius_ then we continue
        if(!circum_sphere_center_opt) continue;
        const vec3 face_b_circum_center = circum_sphere_center_opt.value();

        //Make sure the circumsphere is empty
        if(!circum_sphere_is_empty(face_b_circum_center, point_a, point_b, p)) continue;

        //Compute the angle between the cicumsphere of face_a and the newly created circumsphere
        double theta = ball_pivoting_algorithm::compute_angle_between_circumspheres(
            face_a_circum_center,
            face_b_circum_center,
            edge_center_vec
        );


        //If the node p, creates the most flat surface wrt to face_a, then it is the best 
        //candidate
        if(theta < min_theta){
            min_theta = theta;
            best_candidate = p;
        }
    }

    return best_candidate;
}
//----------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------
//Returns if the point is already part of the triangulated region of the cell surface
bool ball_pivoting_algorithm::node_is_inner_vertex(const oriented_point& p) const noexcept{
    return node_is_inner_vertex(p.id_);
}


bool ball_pivoting_algorithm::node_is_inner_vertex(const unsigned point_id) const noexcept{

    //If the node is not connected to any face
    const auto face_it =  node_face_lst_.find(point_id);

    if(face_it == node_face_lst_.end()) return false;
    if(face_it->second.size() == 0)     return false;

    //If the node is not yet connected to any edge
    const auto edge_it =  node_edge_lst_.find(point_id) ;
    if(edge_it == node_edge_lst_.end()) return false;
    if(edge_it->second.size() == 0)     return false;

    //If it is connected to edges, check that none of these edges are at the triangulation front
    return std::all_of(edge_lst_.begin(), edge_lst_.end(), [](const edge& e) -> bool {return e.is_manifold();});
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Make sure that the normal of the face {f1_a, f1_b, p} points in the same direction as the normal
//of the face {f1_a, f1_b, f1_c}
bool ball_pivoting_algorithm::node_is_compatible_with_edge(
        const oriented_point& p, 
        const oriented_point& f1_a,
        const oriented_point& f1_b,
        const oriented_point& f1_c

    ) noexcept{

    const vec3 v1 = f1_b.position_ - f1_a.position_;
    const vec3 v2 = f1_c.position_ - f1_a.position_;
    const vec3 v3 = p.position_    - f1_a.position_;

    //Compute the normals of the triangles
    const vec3 n1 = v1.cross(v2);
    const vec3 n2 = v3.cross(v1);

    return n1.dot(n2) > 0.;
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Compute the angle between the center of the circumspheres and the center of the edge
//the candidate point, with the smallest angle is the one choosen. Basically the 
//algorithm tries to create a flat triangulated surface. The returned angle should be in the 
//range [0, 2pi]
double ball_pivoting_algorithm::compute_angle_between_circumspheres(
        const vec3& face_a_circum_center,
        const vec3& face_b_circum_center,
        const vec3& edge_center_vec
    ) noexcept{

    const vec3 v1 = (face_a_circum_center - edge_center_vec).normalize();
    const vec3 v2 = (face_b_circum_center - edge_center_vec).normalize();
    double theta = std::acos(v1.dot(v2));
    return theta;
}




//----------------------------------------------------------------------------------------------
//Fill the remaining holes in the surface with triangles
void ball_pivoting_algorithm::fill_surface_holes() noexcept(false){

    //Loop over the edges
    for(edge& e: edge_lst_){

        //If the edge is next to a hole
        if(!e.is_manifold()){


            //-------------------------------------------------------------------------------------------------------
            //Get the nodes connected to the edge
            const unsigned n_a_id = e.n1();
            const unsigned n_b_id = e.n2();

            //Keep track of the nodes of the hole
            std::vector<unsigned> hole_node_lst{n_a_id, n_b_id};

            auto previous_node_id = n_a_id;
            auto current_node_id = n_b_id;

            //std::cout << "Hole found" << std::endl;
            //Loop over the hole, holes with a size greater than 10 nodes won't be filled
            for(unsigned i = 0; i < 11; i++){

                const unsigned next_node_id = get_next_node_of_the_hole(previous_node_id, current_node_id);
                
                //If we are back to the first node, then we have finished the hole
                if(next_node_id == n_a_id) break;
                
                //Add the node to the list of nodes of the hole and continue
                hole_node_lst.push_back(next_node_id);
                previous_node_id = current_node_id;
                current_node_id = next_node_id;

                //If the hole is too big, then throw an exception
                if(i == 10) throw bpa_exception("Hole too big in the surface");
            }
            if(hole_node_lst.size() < 3) throw bpa_exception("Hole too small in the surface");
            //-------------------------------------------------------------------------------------------------------


            //Now that we know the nodes that are part of the hole, we can create the triangles
            //First compute the center of the hole
            vec3 hole_center, hole_center_normal;
    

            for(const unsigned node_id: hole_node_lst){
                hole_center = hole_center + node_lst_[node_id].position_;
                hole_center_normal = hole_center_normal + node_lst_[node_id].normal_;

            }
            hole_center = hole_center / static_cast<double>(hole_node_lst.size());
            hole_center_normal = hole_center_normal / static_cast<double>(hole_node_lst.size());

            //Create a point at the center of the hole, the normal of this point is not important
            const oriented_point hole_center_point(node_lst_.size(), hole_center, hole_center_normal);
            node_lst_.push_back(hole_center_point);
            const unsigned n3_id = node_lst_.size() - 1;

            //Loop over the nodes of the hole
            for(unsigned i = hole_node_lst.size()-1, j = 0 ; j < hole_node_lst.size(); i = j++){
                
                const unsigned n1_id  = hole_node_lst[i];
                const unsigned n2_id  = hole_node_lst[j];
                assert(n1_id != n2_id && n1_id != n3_id && n2_id != n3_id);
                assert(n1_id < node_lst_.size() && n2_id < node_lst_.size() && n3_id < node_lst_.size());

                //Create the triangle, with the correct orientation 
                const unsigned face_id = face_lst_.size();

                //Compute the correct and initial normal of the face
                const vec3 initial_normal = (node_lst_[n2_id].position_ - node_lst_[n1_id].position_).cross(node_lst_[n3_id].position_ - node_lst_[n1_id].position_);
                const vec3 correct_normal = compute_oriented_normal(node_lst_[n1_id], node_lst_[n2_id], node_lst_[n3_id]);
                face& f = initial_normal.dot(correct_normal) > 0 ? face_lst_.emplace_back(n1_id, n2_id, n3_id, face_id) : face_lst_.emplace_back(n1_id, n3_id, n2_id, face_id);

                //Update the node-face list
                node_face_lst_[n1_id].push_back(face_id);
                node_face_lst_[n2_id].push_back(face_id);
                node_face_lst_[n3_id].push_back(face_id);

                //Get the edge between n1 and n2
                auto edge_it = std::find(edge_lst_.begin(), edge_lst_.end(), edge(n1_id, n2_id));
                if(edge_it == edge_lst_.end()) throw bpa_exception("Edge not found");

                //Add the face to the edge
                edge_it->add_face(face_id);

                //Create the 2 other edges of the face
                auto e1_id = get_edge(n1_id, n3_id);
                auto e2_id = get_edge(n2_id, n3_id);
                assert(e1_id < edge_lst_.size() && e2_id < edge_lst_.size());

                edge_lst_[e1_id].add_face(face_id);
                edge_lst_[e2_id].add_face(face_id);

                
            }
            //-------------------------------------------------------------------------------------------------------

        }
    }
}
//----------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------
//Get the next node which is part of the hole
unsigned ball_pivoting_algorithm::get_next_node_of_the_hole(const unsigned previous_node, const unsigned current_node) const noexcept(false){
    assert(current_node < node_lst_.size());

    
    //Loop over the edges to which the node is connected
    for(const unsigned edge_id: node_edge_lst_.at(current_node)){
        const edge& e = edge_lst_[edge_id];

        //If the edge is not part of the triangulated surface, then it is part of the hole
        if(!e.is_manifold()){

            //Make sure that the edge is not the previous edge of the hole
            if(e != edge(previous_node, current_node)){
                return e.n1() == current_node ? e.n2() : e.n1();
            }
        }
    }

    throw bpa_exception("The ball pivoting algortihm has failed to generate a watertight surface");
}
//----------------------------------------------------------------------------------------------


