cd NM_COURSE_WORK
cd c_lib

dir=cmake_build
mkdir $dir
cp main.cpp $dir/
cp CMakeLists.txt $dir/
cd $dir

cmake CMakeLists.txt
make

cd ..
cp cmake_build/nm_course_work nm_course_work