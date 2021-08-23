set(FILES
	"Augmenter.h"
	"Lidar.h"
	"Lisa.h"
	"Material.h"
	"MiniLisa.h"
	"ParticleDist.h"
	"Utils.h"
)
set(RELATIVE_PATH "include")

# Add files
if(FILES)
	#if files is not empty:
	target_sources(pylisa PRIVATE ${FILES})
	
	# if we want to install files
	# this will be installed to CMAKE_INSTALL_PREFIX/include
	install(FILES ${FILES} DESTINATION ${RELATIVE_PATH})
endif()