from spack.package import *


class Amp(CMakePackage):

    homepage = "https://re-git.lanl.gov/xcap/oss/solvers/amp"
    git = "ssh://git@re-git.lanl.gov:10022/xcap/oss/solvers/amp.git"
    
    version("master", branch="master")
    version("3.0.0", tag="3.0.0")

    
    variant("mpi", default=True, description="build with mpi")
    variant("stacktrace", default=False)
    variant("hypre", default=False)
    variant("cuda", default=False)
    variant("rocm", default=False)


    depends_on("cmake@3.26.0:", type="build")
    depends_on("mpi", when="+mpi")
    depends_on("tpl-builder@master")

    tpl_depends = ["stacktrace", "hypre", "cuda", "rocm"]

    for _variant in tpl_depends:
        depends_on("tpl-builder@master+" + _variant, when="+" + _variant)
        depends_on("tpl-builder@master~" + _variant, when="~" + _variant)

    def cmake_args(self):

        args = [
            "-D TPL_DIRECTORY="+self.spec["tpl-builder"].prefix,
            "-D USE_CUDA=0",
            "-D AMP_INSTALL_DIR="+self.spec.prefix,
            "-D CXX_STD=17",
            "-D DISABLE_ALL_TESTS=ON",
        ]
        
        #TODO have amp_data as a dependencie
        

        if self.spec.satisfies("^mpi"):
            args.append("-DCMAKE_CXX_COMPILER=mpicxx")
        return args


