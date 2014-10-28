#ifndef NODE_H
#define NODE_H

#include <utility>
#include <string>
#include <vector>

class Node : public std::pair<double, double>
{
public:
        Node(double x, double y, unsigned int index) : std::pair<double, double>(y,x), index(index)
        {}

        double GetLongitude() const {return second;}
        double GetLatitude() const {return first;}
        
        unsigned int index;
        
        std::vector<std::string> tags;
        std::string tag;
};

#endif /* NODE_H */
