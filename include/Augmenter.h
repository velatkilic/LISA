#ifndef AUGMENTER_H
#define AUGMENTER_H

#include <vector>

class Augmenter {
public:
	//! @brief Abstract augmenter class. Any augmenter should be able to take a point cloud scan and return the augmented scan
	//! @param pc : Point cloud
	//! @return Augmented point cloud
	virtual std::vector<std::vector<double>> augment(std::vector<std::vector<double>>& pc) = 0;
};

#endif // AUGMENTER_H