#include "KDTree.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstring>

KDNode::KDNode(double la, double lo, const std::string &desc) {
    left = NULL;
    right = NULL;
    name = desc;
    lat = la;
    lon = lo;
}

KDNode::~KDNode() {
}

/* provided by instructor */
double KDNode::distance(double _la, double _lo) {
    double param = M_PI / 180.0; // required for conversion from degrees to radians
    double rad = 3956.0;  // radius of earth in miles
    double d_lat = (_la - lat) * param;
    double d_lon = (_lo - lon) * param;
    double dist = sin(d_lat/2) * sin(d_lat/2) + cos(lat*param) * cos(_la*param) * sin(d_lon/2) * sin(d_lon/2);
    dist = 2.0 * atan2(sqrt(dist), sqrt(1-dist));
    return rad * dist;
}

KDTree::KDTree() {
    size = 0;
    root = NULL;
}

KDTree::~KDTree() {
    destroy(root);
}

void KDTree::destroy(KDNode *p) {
    if (p) {    // postorder traversal
        destroy(p->left);
        destroy(p->right);
        delete p;
    }
}

// return the number of nodes in the tree
unsigned int KDTree::getSize() {
    return size;
}

// insert a new node into the tree
void KDTree::insert(double la, double lo, const std::string &desc) {
    // pointer to root of tree with an initial depth of 0
    root = insert(root, 0, la, lo, desc);
}

// recursively insert a new node k-dimensionally (partitioning) into the tree
KDNode *KDTree::insert(KDNode *p, int depth, double la, double lo, const std::string &desc) {
    if (!p) {                               // node is a null
        size++;                             // increment number of nodes in the tree
        return new KDNode(la, lo, desc);    // return a new node
    }
    
    // node has a depth of 0 -- latitude partitioning
    // node has a depth of 1 -- longitude partitioning
    switch (depth) {
        case 0:
            if (la >= p->lat) {     // new latitude is greater than or equal to node->latitude
                p->right = insert(p->right, 1, la, lo, desc);
            } else {                // new latitude is less than node->latitude
                p->left = insert(p->left, 1, la, lo, desc);
            }
            break;
        case 1:
            if (lo >= p->lon) {     // new longitude is greater than or equal to node->longitude
                p->right = insert(p->right, 0, la, lo, desc);
            } else {                // new longitude is less than node->longitude
                p->left = insert(p->left, 0, la, lo, desc);
            }
            break;
    }
    return p;   // return node
}

// print all the nodes under a distance 'rad' from 'la,lo' and where
// filter is a non-empty substring of their description
// data type changed to return a vector of kdnode pointers
std::vector <KDNode*> KDTree::printNeighbors(double la, double lo, double rad, const std::string &filter) {
    
    //using a stack and vector to store data structure
    std::stack<KDNode*> query;      // holds nodes as they are found they are found in rangeQuery
    std::vector<KDNode*> nodes;
    std::ofstream file;
    unsigned int count = 0, comparisons = 0;
    
    rangeQuery(root, &query, &count, &comparisons, 0, la, lo, rad, filter);
    
    //open the file that contains the points
    file.open("data/markers.js");
    
    file << "var markers = [\n";
    file << "\t[\"" << "CENTER" << "\", " << la << ", " << lo << "],\n";
    // from the stack, move it to a vector
    nodes = std::vector<KDNode*>(count);
    for (int i = 0; i < count; i++) {
    	nodes[i] = query.top();
    	file << "\t[\"" << nodes[i]->name << "\", " << nodes[i]->lat << ", " << nodes[i]->lon << "],\n";
    	query.pop();
    }
    file << "];\n";
    file.close();
    
    std::cerr << comparisons << " comparisons were made\n";
    
    return nodes;
}

// recursively prune and print a node under a distance 'rad' from 'la, lo'
void KDTree::rangeQuery(KDNode *p, std::stack<KDNode*> *query, unsigned int *count, unsigned int *comparisons, int depth, double la, double lo, double rad, const std::string &filter) {
    if (!p) {   // node is null
        return;
    }
    // if the distance between the node and the center is less than the radius, it lies within the radius
    else if (p->distance(la,lo) < rad && p->name.find(filter) != std::string::npos) {
        
        //pushes the points into a stack
        query->push(p);
        *count += 1;    // increment number of records fetched
    }
    
    *comparisons += 1;  // increment number of comparisions made
    
    // node has a depth of 0 -- latitude partitioning
    // node has a depth of 1 -- longitude partitioning
    switch (depth) {    // ...pruning...
        case 0:
            if (p->distance(la, p->lon) > rad && la < p->lat)           // if farthest radial point east of center is less than node->latitude
                rangeQuery(p->left, query, count, comparisons, 1, la, lo, rad, filter);
            else if (p->distance(la, p->lon) >= rad && la > p->lat)     // if farthest radial point west of center is greater than or equal to node->latitude
                rangeQuery(p->right, query, count, comparisons, 1, la, lo, rad, filter);
            else {                                                      // circular perimeter lies on both sides of partition
                rangeQuery(p->left, query, count, comparisons, 1, la, lo, rad, filter);
                rangeQuery(p->right, query, count, comparisons, 1, la, lo, rad, filter);
            }
            break;
        case 1:
            if (p->distance(p->lat, lo) > rad && lo < p->lon)           // if farthest radial point south of center is less than node->longitude
                rangeQuery(p->left, query, count, comparisons, 0, la, lo, rad, filter);
            else if (p->distance(p->lat, lo) >= rad && lo > p->lon)     // if farthest radial point north of center is greater than or equal to node->longitude
                rangeQuery(p->right, query, count, comparisons, 0, la, lo, rad, filter);
            else {                                                      // circular perimeter lies on both sides of partition
                rangeQuery(p->left, query, count, comparisons, 0, la, lo, rad, filter);
                rangeQuery(p->right, query, count, comparisons, 0, la, lo, rad, filter);
            }
            break;
    }
}

// merge sort
void KDTree::merge(std::vector<KDNode*>& A, std::vector<KDNode*>& aux, ul_int lo, ul_int mid, ul_int hi) {
    // copy array
    std::memcpy(&aux[0] + lo, &A[0] + lo, (hi-lo+1)*sizeof(KDNode*));

    // merge
    ul_int i = lo, j = mid + 1;
    for (int k = lo; k <= hi; k++) {
        if (i > mid) A[k] = aux[j++];
        else if (j > hi) A[k] = aux[i++];
        else if (   // checks if point aux[j] has a smaller angle to A[0] than aux[i] to A[0]
        	atan2(aux[j]->lat - A[0]->lat, aux[j]->lon - A[0]->lon) * 180 / M_PI <
        	atan2(aux[i]->lat - A[0]->lat, aux[i]->lon - A[0]->lon) * 180 / M_PI
        	) A[k] = aux[j++];
        else A[k] = aux[i++];
    }
}

void KDTree::r_mergesort(std::vector<KDNode*>& A, std::vector<KDNode*>& aux, ul_int lo, ul_int hi) {
    // base case
    if (hi <= lo) return;
    // divide
    ul_int mid = lo + (hi - lo) / 2;
    // recursively sort halves
    r_mergesort(A, aux, lo, mid);
    r_mergesort(A, aux, mid+1, hi);
    // merge results
    merge(A, aux, lo, mid, hi);
}

void KDTree::mergesort(std::vector<KDNode*>& A, ul_int n) {
    // allocate space for aux
    std::vector<KDNode*> aux (n);
    // call recursive mergesort
    r_mergesort(A, aux, 1, n - 1);  // begin at first index
}

// checks the orientation of the ponts to make sure they make a right turn
double KDTree::cwturn(KDNode *a, KDNode *b, KDNode *c){
    double val = ((b->lat-a->lat)*(c->lon-b->lon))-((b->lon - a->lon)*(c->lat - b->lat));
    if (val==0) return 0;   // colinear 
    return(val>0)?1:2;      // the turn
}

// grabs the node that is right under the top
KDNode * KDTree::nextToTop(std::stack<KDNode*> &v_stack){
    KDNode* p = v_stack.top();
	v_stack.pop();
	KDNode* res = v_stack.top();
	v_stack.push(p);
	return res;
}

// outputs convex points of a set of points passed from printNeighbor
unsigned int KDTree::printConvexHull(double la, double lo, double rad, const std::string &filter) {
    // create vector to hold the nodes
    std::vector<KDNode*> nodes;
    std::ofstream file;
    // printNeighbor will return a vector of successful results
    nodes = printNeighbors(la, lo, rad, filter);
    
    KDNode *temp;
    int i, j;
    // insertion sort -- from lowest lattitude
    for (i = 0; i < nodes.size(); i++) {
        temp = nodes[i];
        // inserts A[j] into the right place in sorted part
        for (j = i; j > 0 && nodes[j-1]->lat > temp->lat; j--) {
            nodes[j] = nodes[j-1];
        }
        nodes[j] = temp;
    }
    
    // sort by polar angles
    mergesort(nodes, nodes.size());
    
	// print angles to make sure in-order
    for (int i = 0; i < nodes.size(); i++)
    	std::cerr << "\t[\"" << nodes[i]->name << "\", " << nodes[i]->lat << ", " << nodes[i]->lon << "],\t" << atan2(nodes[i]->lat - nodes[0]->lat, nodes[i]->lon - nodes[0]->lon) * 180 / M_PI << "\n";
    
    //sets the anchor, lowest y value
    KDNode *anchor = nodes[0];
    int n = 1;
    
    //graham scan
    for(int i = 1 ; i < nodes.size(); i++) {
        while(i < nodes.size() - 1 && cwturn(anchor, nodes[i], nodes[i + 1]) == 0) {
            i++;
        }
        nodes[n] = nodes[i];
        n++;
    }
    
    //convex hull cant be made if size is less than 3
    if(n < 3) return nodes.size();

    //pushes the ordered points to a stack
    std::stack<KDNode*> v_stack;
    v_stack.push(nodes[0]);
    v_stack.push(nodes[1]);
    v_stack.push(nodes[2]);
    
    //passes the ordered points to a helper function
    //deletes a node if the turn isn't right
    for (int i = 3; i < n; i++) {
	    while (cwturn(nextToTop(v_stack), v_stack.top(), nodes[i]) != 2)
		    v_stack.pop();
	    v_stack.push(nodes[i]);
    }
    
    file.open("data/convex.js");
    
    file << "var convex = [\n"  ;
    //prints the convex points
    while(!v_stack.empty()) {
        KDNode *v = v_stack.top();
        file << "\t{lat: " << v->lat << ", lng: " << v->lon <<"},\n";
        v_stack.pop();
    }
    file << "];\n";
    file.close();
    
    return nodes.size();
}

