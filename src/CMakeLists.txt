set(FILES 
	"Lidar.cpp"
	"Lisa.cpp"
	"MiniLisa.cpp"
	"ParticleDist.cpp"
	"Utils.cpp"
)
set(RELATIVE_PATH "src")

# Add files
if(FILES)
	#if files is not empty:
	target_sources(pylisa PRIVATE ${FILES})
	
	# if we want to install files:
	# this will be installed to CMAKE_INSTALL_PREFIX/src
	install(FILES ${FILES} DESTINATION ${RELATIVE_PATH})
endif()