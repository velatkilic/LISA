#include"Lidar.h"

//! @brief Overload << operator to print lidar object fields to standard out
//! @param os : Output stream
//! @param lidar : Lidar object
//! @return Output stream object
std::ostream& operator<<(std::ostream& os, const Lidar& lidar) {
	os << "Range uncertainty [m] " << lidar.ran_uncer << std::endl
		<< "Max range [m] " << lidar.max_range << std::endl
		<< "Min range [m] " << lidar.min_range << std::endl
		<< "Laser wavelength [um] " << lidar.wavelength << std::endl;
	return os;
}