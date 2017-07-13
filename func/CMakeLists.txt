#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( BASENAME_ONLY
		  LIBRARY_NAME dune_GammaTagging
		  LIB_LIBRARIES larcore_Geometry_Geometry_service
                        larcore_Geometry
                        larData_Utilities
						${ART_FRAMEWORK_SERVICE_REGISTRY}
						${CETLIB}
						${ROOT_BASIC_LIB_LIST}
						${Boost_FILESYSTEM_LIBRARY}
						${Boost_SYSTEM_LIBRARY}
		 DICT_LIBRARIES lardataobj_RecoBase
						

)

install_headers()
install_source()
install_fhicl()
