from spack.package import *


class Amp(CMakePackage):

    homepage = "https://re-git.lanl.gov/xcap/oss/solvers/amp"
    git = "ssh://git@re-git.lanl.gov:10022/xcap/oss/solvers/amp.git"
    
    version("master", branch="master")
    version("3.0.0", tag="3.0.0")

    
    variant("mpi", default=True, description="build with mpi")
    variant("hypre", default=False)
    variant("cuda", default=False)
    variant("rocm", default=False)
    variant("openmp", default=False)

    depends_on("stacktrace@master")
    depends_on("cmake@3.26.0:", type="build")
    depends_on("mpi", when="+mpi")
    depends_on("tpl-builder@master+stacktrace")

    conflicts("+rocm +cuda")

    tpl_depends = ["hypre", "cuda", "rocm"]


    for _variant in tpl_depends:
        depends_on("tpl-builder@master+" + _variant, when="+" + _variant)
        depends_on("tpl-builder@master~" + _variant, when="~" + _variant)



    def cmake_args(self):

        args = [
            "-D TPL_DIRECTORY="+self.spec["tpl-builder"].prefix,
            "-D USE_HIP={condition}".format(condition="1" if self.spec.satisfies("+rocm") else "0"),         
            "-D AMP_INSTALL_DIR="+self.spec.prefix,
            "-D CXX_STD=17",
            "-D DISABLE_ALL_TESTS=ON",
            self.define_from_variant("USE_OPENMP", "openmp"),
            self.degine_from_variant("USE_CUDA", "cuda")
        ]
        
        #TODO have amp_data as a dependencie
        
        if self.spec.satisfies("+mpi"):
            args.append("-DCMAKE_CXX_COMPILER=mpicxx")
        
        if self.spec.satisfies("+cuda"):
            args.append("-DCMAKE_CUDA_FLAGS=--extended-lambda") 

 
        return args
