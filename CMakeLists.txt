simple_plugin(gammatagging "module"
                        larcore_Geometry_Geometry_service
                        larcore_Geometry
                        larData_Utilities
                        lardataobj_RecoBase
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        
                        nusimdata_SimulationBase
                        larsim_MCCheater_PhotonBackTracker_service
						
                        lardataobj_RawData
						larsim_MCCheater_BackTracker_service
					
						${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
						${ART_FRAMEWORK_PRINCIPAL}
						${ART_FRAMEWORK_SERVICES_REGISTRY}
						${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
						art_Persistency_Common canvas
						art_Persistency_Provenance canvas
						art_Utilities canvas
						${MF_MESSAGELOGGER}
						${MF_UTILITIES}
						cetlib cetlib_except
						${ROOT_BASIC_LIB_LIST}
				BASENAME_ONLY
)

add_subdirectory(func)

install_headers()
install_fhicl()
install_source()
install_scripts()
