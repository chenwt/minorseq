# global properties
GET_PROPERTY(LOCAL_COMPILE_FLAGS GLOBAL PROPERTY MINORSEQ_COMPILE_FLAGS_GLOBAL)
GET_PROPERTY(LOCAL_LINK_FLAGS GLOBAL PROPERTY MINORSEQ_LINK_FLAGS_GLOBAL)

# get all header files for IDE support
file(GLOB_RECURSE MS_HEADER "${MS_IncludeDir}/*.h")
file(GLOB_RECURSE MS_HIDDEN_HEADER "*.h")
# get sources for src/
file(GLOB MS_CPP        "*.cpp")
# get sources for src/tools
file(GLOB MS_TOOLS_CPP  "tools/*.cpp")

# ssw C library cannot be build with ${LOCAL_COMPILE_FLAGS} since it's C
add_library(ssw STATIC ${ssw_INCLUDE_DIRS}/ssw.c)

# add main library including everything
add_library(minorseq STATIC
    ${MS_HEADER}
    ${MS_HIDDEN_HEADER}
    ${ssw_INCLUDE_DIRS}/ssw_cpp.cpp
    ${MS_CPP}
)

target_link_libraries(minorseq
    ${ZLIB_LDFLAGS}
    ${HTSLIB_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    ${pbcopper_LIBRARIES}
    ssw
)

target_include_directories(minorseq PUBLIC
    ${MS_IncludeDir}
    ${CMAKE_BINARY_DIR}/generated
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${PacBioBAM_INCLUDE_DIRS}
    ${pbcopper_INCLUDE_DIRS}
    ${ssw_INCLUDE_DIRS}
)

set_target_properties(minorseq PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILE_FLAGS})
if (LOCAL_LINK_FLAGS)
    set_target_properties(minorseq PROPERTIES LINK_FLAGS ${LOCAL_LINK_FLAGS})
endif()

if (MS_build_bin)
    function(create_exe exeName)
        add_executable(${exeName} ${MS_SourceDir}/tools/main/${exeName}.cpp)

        set_target_properties(${exeName} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        )

        target_link_libraries(${exeName}
            ${CMAKE_THREAD_LIBS_INIT}
            ${CMAKE_DL_LIBS}
            minorseqtools
        )

        install(TARGETS ${exeName} RUNTIME DESTINATION bin)
        set_target_properties(${exeName} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILE_FLAGS})
        if (LOCAL_LINK_FLAGS)
            set_target_properties(${exeName} PROPERTIES LINK_FLAGS ${LOCAL_LINK_FLAGS})
        endif()
    endfunction()

    add_library(minorseqtools STATIC
        ${MS_TOOLS_CPP}
    )

    target_link_libraries(minorseqtools minorseq)

    set_target_properties(minorseqtools PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILE_FLAGS})
    if (LOCAL_LINK_FLAGS)
        set_target_properties(minorseqtools PROPERTIES LINK_FLAGS ${LOCAL_LINK_FLAGS})
    endif()

    create_exe(juliet)
    create_exe(fuse)
    create_exe(cleric)
endif()
