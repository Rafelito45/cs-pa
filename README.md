# cs-pa

### Assignment

Final programming assignment for Data Structures and Algorithms course

http://homepage.cs.uri.edu/~malvarez/teaching/2016s/csc212/pa-3.html

Complete source files provided:
* **data/index.html**
* **data/markers.js**
* **data/convex.js**
* **data/parse-osm.py**
* **data/boston_massachusetts.txt**
* **search-map.cpp**  -- *slightly modified and renamed* **search-map-kd.cpp**

### Group members
* Susallin Chean
* Dion Rubio
* Osto Vargas (***me!***)

Although this was a group effort, I have contributed a substantial amount of code and logic into the **KDTree.h** and **KDTree.cpp** source files.

### My contributions

* Project lead
* Entire **KDTree.h** source file -- wrote all necessary methods data members and expected outcomes
* Post order deletion of **KDTree**
* Merge sort implementation for sorting an array of **KDNodes** by polar angle to a fixed hull -- ***O(nlogn)***
* Insertion of **KDNode** objects within a **KDTree**
* Switch statements utilizing depth of **KDTree** traversal for partition space pruning in ***KDTree::insert*** & ***KDTree::rangeQuery***
* ***std::stack<T>*** implementation to store matching results on the fly as they are found in **KDTree::printNeighbors**
* ***std::ofstream*** output of **markers.js** and **convex.js** files
* Pointer guru of the group -- handled passing of vectors, stacks, and other data types within method parameters
* Documentation and *Refractoring of source files

*After the assignment was handed in I refractored the logic behind the convex hull's counter clockwise algorithm written by Dion Rubio on my free time. I found a great source example written in Java by Robert Sedgewick and Kevin Wayne.
