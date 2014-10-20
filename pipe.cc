#include "pipe.h"

World worldFromString(std::string data) {
	
	std::stringstream d(data);
	std::string line, entry, name, x, y, z;
	
	while (std::getline(d, line))
	{
		std::stringstream l(line);
		std::getline(l, entry, ' ');

		if (entry.compare("s") == 0)
		{
			std::getline(l, name, ' ');
			std::getline(l, x, ' ');
			std::getline(l, y, ' ');
			std::getline(l, z, ' ');
		}
		else if (entry.compare("c") == 0)
		{
			getline(l, name, ' ');
		}
		else if (entry.compare("r") == 0)
		{
			getline(l, name, ' ');
		}
	}
}
