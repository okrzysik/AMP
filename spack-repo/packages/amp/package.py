from spack.package import *


class Amp(CMakePackage):

    homepage = "https://re-git.lanl.gov/xcap/oss/solvers/amp"
    git = "ssh://git@re-git.lanl.gov:10022/xcap/oss/solvers/amp.git"
    
    version("master", branch="master")
    version("3.0.0", tag="3.0.0")

    
    variant("mpi", default=True, description="build with mpi")
    variant("hypre", default=False)
    variant("cuda", default=False)
    variant("cuda_arch", default="none", values = ("none", "10", "11", "12", "13", "20", "21", "30", "32", "35", "37", "50", "52", "53", "60", "61", "62", "70", "72", "75", "80", "86", "87", "89", "90"), multi=False)
    variant("rocm", default=False)

    depends_on("stacktrace@master")
    depends_on("cmake@3.26.0:", type="build")
    depends_on("mpi", when="+mpi")
    depends_on("tpl-builder@master+stacktrace")

    conflicts("+rocm +cuda")

    tpl_depends = ["hypre", "cuda", "rocm"]


    for _variant in tpl_depends:
        depends_on("tpl-builder@master+" + _variant, when="+" + _variant)
        depends_on("tpl-builder@master~" + _variant, when="~" + _variant)

    #shamelessly stollen from the hypre spack package
    for sm_ in CudaPackage.cuda_arch_values:
        depends_on(
            "tpl-builder@master+cuda cuda_arch={0}".format(sm_),
            when="+cuda cuda_arch={0}".format(sm_),
        )


    def cmake_args(self):

        args = [
            "-D TPL_DIRECTORY="+self.spec["tpl-builder"].prefix,
            "-D USE_HIP={condition}".format(condition="1" if self.spec.satisfies("+rocm") else "0"),         
            "-D AMP_INSTALL_DIR="+self.spec.prefix,
            "-D CXX_STD=17",
            "-D DISABLE_ALL_TESTS=ON",
        ]
        
        #TODO have amp_data as a dependencie
        
        if self.spec.satisfies("+mpi"):
            args.append("-DCMAKE_CXX_COMPILER=mpicxx")
        return args
        
        if self.spec.satisfies("+cuda"):
            args.append("-D USE_CUDA=1")
            args.append=("-D CMAKE_CUDA_ARCHITECTURES=70")
            #args.append("-D CUDA_ARCHITECTURES=" + self.spec.variants["cuda_arch"].value)
            args.append("-D CMAKE_CUDA_COMPILER=" + self.spec["cuda"].prefix + "bin/nvcc")

        else:
            args.append("-D USE_CUDA=0")
