CXX      = g++
CXXFLAGS = -O3 -Wno-deprecated
#CXXFLAGS = -O3 -xP -Wno-deprecated

all: VolumeRendering

LevelSet3DView:
	${CXX} ${CXXFLAGS} ${VTK_INCLUDE} LevelSet3DView.cxx ${VTK_LIB}
	./a.out /data/sims/sim3D_8/mycell.19999.Real64

VolumeRendering:
	${CXX} ${CXXFLAGS} ${VTK_INCLUDE} VolumeRendering.cxx ${VTK_LIB} 
	./a.out "/tmp/phi.214.Real64" \
					"/tmp/phi.214.jpg"

VolumeRenderingGrayScott.o:
	${CXX} ${CXXFLAGS} ${VTK_INCLUDE} VolumeRenderingGrayScott.cxx ${VTK_LIB} 
	./a.out "/home/abergman/Research/DCell/temp/vec.1190.dat" \
					"/home/abergman/Research/DCell/temp/vec.1190.jpg"

TransparentContourGrayScott.o:
	${CXX} ${CXXFLAGS} ${VTK_INCLUDE} TransparentContourGrayScott.cxx ${VTK_LIB} 
	./a.out "/home/abergman/Research/DCell/temp/vec.1190.dat" \
					"/home/abergman/Research/DCell/temp/vec.1190.jpg"

OffScreenCone.o:
	${CXX} ${CXXFLAGS} ${VTK_INCLUDE} OffScreenCone.cxx ${VTK_LIB}
	./a.out 

sync:
	rsync -av --delete /home/abergman/Research/DCell/VTK/ levlabhn:/home/abergman/Research/DCell/VTK/
	ssh 
	

submit:
	scp a.out levlabhn:exec/
	qsub.sh -t 1-333 ~/Research/sge/BatchProcess.sh a.out /home/share/sim chem.III.Real64 output.III.jpg
	mencoder "mf://*.?.jpg" "mf://*.??.jpg" "mf://*.???.jpg" -quiet -mf fps=30 -o chem.avi -ovc lavc -lavcopts vcodec=msmpeg4v2
	mencoder "mf://*.??.jpg" "mf://*.???.jpg" "mf://*.????.jpg" "mf://*.?????.jpg" -quiet -mf fps=30 -o chem.avi -ovc lavc -lavcopts vcodec=wmv2:vbitrate=2000

ifeq (${HOSTNAME},"levlabhn.clark.jhu.edu")
#OFF-SCREEN LIBRARIES
  VTK_DIR=/share/apps/VTK
  VTK_INCLUDE= \
-I${VTK_DIR} \
-I${VTK_DIR}/Utilities \
-I${VTK_DIR}/VolumeRendering \
-I${VTK_DIR}/Rendering \
-I${VTK_DIR}/Hybrid \
-I${VTK_DIR}/Widgets \
-I${VTK_DIR}/Rendering/Testing/Cxx \
-I${VTK_DIR}/IO \
-I${VTK_DIR}/Imaging \
-I${VTK_DIR}/Graphics \
-I${VTK_DIR}/GenericFiltering \
-I${VTK_DIR}/Filtering \
-I${VTK_DIR}/Common \
-I${VTK_DIR}/Common/Testing/Cxx \
-I${VTK_DIR}/Utilities/DICOMParser \
-I${VTK_DIR}/Utilities/vtkfreetype/include \
-I${VTK_DIR}/Utilities/vtknetcdf \
-I${VTK_DIR}/Utilities/vtkexodus2/include
  VTK_LIB=-L${VTK_DIR}/bin -Wl,-rpath,${VTK_DIR}/bin \
-lvtkHybrid   -lvtkmetaio      -lvtkRendering        -lvtkWidgets \
-lmpistubs        -lvtkFiltering         -lvtkImaging  -lvtksqlite           -lvtkzlib \
-lvtkalglib       -lvtkfreetype          -lvtkInfovis  -lvtksys \
-lvtkCharts       -lvtkftgl              -lvtkIO       -lvtktiff \
-lvtkCommon       -lvtkGenericFiltering  -lvtkjpeg     -lvtkNetCDF      -lvtkverdict \
-lvtkDICOMParser  -lvtklibxml2  -lvtkpng         -lvtkViews \
-lvtkexoIIc       -lvtkGraphics          -lvtkproj4       -lvtkVolumeRendering
else
#ON-SCREEN LIBRARIES
VTK_DIR=/share/apps/vtk/VTK
VTK_INCLUDE= \
-I${VTK_DIR} \
-I${VTK_DIR}/Utilities \
-I${VTK_DIR}/VolumeRendering \
-I${VTK_DIR}/Rendering \
-I${VTK_DIR}/Hybrid \
-I${VTK_DIR}/Widgets \
-I${VTK_DIR}/Rendering/Testing/Cxx \
-I${VTK_DIR}/IO \
-I${VTK_DIR}/Imaging \
-I${VTK_DIR}/Graphics \
-I${VTK_DIR}/GenericFiltering \
-I${VTK_DIR}/Filtering \
-I${VTK_DIR}/Common \
-I${VTK_DIR}/Common/Testing/Cxx \
-I${VTK_DIR}/Utilities/DICOMParser \
-I${VTK_DIR}/Utilities/vtkfreetype/include \
-I${VTK_DIR}/Utilities/vtknetcdf \
-I${VTK_DIR}/Utilities/vtkexodus2/include
VTK_LIB=-L${VTK_DIR}/bin -Wl,-rpath,${VTK_DIR}/bin \
-lvtkHybrid       -lvtkmetaio            -lvtkRendering  -lvtkWidgets \
-lmpistubs        -lvtkFiltering         -lvtkImaging    -lvtksqlite           -lvtkzlib \
-lvtkalglib       -lvtkfreetype          -lvtkInfovis    -lvtksys \
-lvtkCharts       -lvtkftgl              -lvtkIO         -lvtktiff \
-lvtkCommon       -lvtkGenericFiltering  -lvtkjpeg       -lvtkNetCDF      -lvtkverdict \
-lvtkDICOMParser  -lvtklibxml2           -lvtkpng        -lvtkViews \
-lvtkexoIIc       -lvtkGraphics          -lvtkproj4      -lvtkVolumeRendering
endif

OBJS = LevelSetView.o

Medical1.o: 
	icpc -Wno-deprecated -DNDEBUG -c Medical1.cxx ${VTK_INCLUDE}

LevelSetView.o: LevelSetView.cxx
	icpc -O3 -xP -Wno-deprecated -c LevelSetView.cxx ${VTK_INCLUDE}

link: ${OBJS}
	icpc -O3 -xP  -no-ipo ${OBJS} ${VTK_LIB}

linkMedical1: Medical1.o
	icpc -Wno-deprecated -O3 -DNDEBUG  -fPIC -xP -no-ipo Medical1.o \
	-o Medical1 -rdynamic \
	-L/usr/lib/gcc/x86_64-redhat-linux/3.4.3 \
	${VTK_LIB}

runMedical1: linkMedical1
	./Medical1  ~/ProgramFiles/VTKData/Data/headsq/quarter

run: link
	./a.out

.PHONY: run link all
