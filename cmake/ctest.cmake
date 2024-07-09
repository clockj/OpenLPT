# Enable testing with CTest
enable_testing()

# test header
set(TESTHEADERFILES
    ${CMAKE_HOME_DIRECTORY}/test/test.h
)


# test Matrix
add_executable(test_Matrix 
    ${CMAKE_HOME_DIRECTORY}/test/test_Matrix.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Matrix PRIVATE Matrix)
add_test(test_Matrix test_Matrix)
set_tests_properties(test_Matrix
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test myMATH
add_executable(test_myMATH 
    ${CMAKE_HOME_DIRECTORY}/test/test_myMATH.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_myMATH PRIVATE myMath Matrix)
add_test(test_myMATH test_myMATH)
set_tests_properties(test_myMATH
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ImageIO
add_executable(test_ImageIO 
    ${CMAKE_HOME_DIRECTORY}/test/test_ImageIO.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ImageIO PRIVATE ImageIO Matrix tiff)
add_test(test_ImageIO test_ImageIO)
set_tests_properties(test_ImageIO
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test Camera
add_executable(test_Camera 
    ${CMAKE_HOME_DIRECTORY}/test/test_Camera.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Camera PRIVATE Camera Matrix myMath)
add_test(test_Camera test_Camera)
set_tests_properties(test_Camera
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectInfo
add_executable(test_ObjectInfo 
    ${CMAKE_HOME_DIRECTORY}/test/test_ObjectInfo.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ObjectInfo PRIVATE ObjectInfo Camera)
add_test(test_ObjectInfo test_ObjectInfo)
set_tests_properties(test_ObjectInfo
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")
# test SphereInfo
add_executable(test_SphereInfo 
    ${CMAKE_HOME_DIRECTORY}/test/test_SphereInfo.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_SphereInfo PRIVATE SphereInfo Camera)
add_test(test_SphereInfo test_SphereInfo)
set_tests_properties(test_SphereInfo
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectFinder
add_executable(test_ObjectFinder 
    ${CMAKE_HOME_DIRECTORY}/test/test_ObjectFinder.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ObjectFinder PRIVATE ObjectFinder ObjectInfo ImageIO)
add_test(test_ObjectFinder test_ObjectFinder)
set_tests_properties(test_ObjectFinder
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test StereoMatch
add_executable(test_StereoMatch 
    ${CMAKE_HOME_DIRECTORY}/test/test_StereoMatch.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_StereoMatch PRIVATE Camera StereoMatch ObjectInfo ImageIO ObjectFinder myMath Matrix)
add_test(test_StereoMatch test_StereoMatch)
set_tests_properties(test_StereoMatch
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test OTF
add_executable(test_OTF 
    ${CMAKE_HOME_DIRECTORY}/test/test_OTF.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_OTF PRIVATE OTF Matrix)
add_test(test_OTF test_OTF)
set_tests_properties(test_OTF
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test Shake
add_executable(test_Shake 
    ${CMAKE_HOME_DIRECTORY}/test/test_Shake.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Shake PRIVATE Shake Shake ObjectFinder StereoMatch ObjectInfo Camera ImageIO OTF myMath Matrix)
add_test(test_Shake test_Shake)
set_tests_properties(test_Shake
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test IPR
add_executable(test_IPR 
    ${CMAKE_HOME_DIRECTORY}/test/test_IPR.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_IPR PRIVATE 
    Matrix
    myMath
    Camera
    ImageIO
    ObjectInfo
    ObjectFinder
    StereoMatch
    Shake
    OTF
    IPR
)
add_test(test_IPR test_IPR)
set_tests_properties(test_IPR
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test PredField
add_executable(test_PredField 
    ${CMAKE_HOME_DIRECTORY}/test/test_PredField.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_PredField PRIVATE PredField Matrix ObjectInfo)
add_test(test_PredField test_PredField)


# test Track
add_executable(test_Track 
    ${CMAKE_HOME_DIRECTORY}/test/test_Track.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Track PRIVATE Track Matrix ObjectInfo)
add_test(test_Track test_Track)


# test STB
add_executable(test_STB 
    ${CMAKE_HOME_DIRECTORY}/test/test_STB.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_STB PRIVATE ImageIO STB Matrix ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track)
add_test(test_STB test_STB)


# test PIVChallenge
add_executable(test_PIVChallenge
    ${CMAKE_HOME_DIRECTORY}/test/test_PIVChallenge.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_PIVChallenge PRIVATE ImageIO STB Matrix ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track)
add_test(test_PIVChallenge test_PIVChallenge)