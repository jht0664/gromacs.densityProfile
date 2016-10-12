# gromacs.util.densityProfile
## About
Fortran code for Density Profile program with Gromacs 4.x and 5.x

## Requirement
  * installed xdr library supported by Frans van Hoesel.
    * See [GibHub Page](https://github.com/Pappulab/xdrf) to get detail.
    * For convenience, you can get the v1.2 'libxdrf_v1.2.tgz', in **xdr** folder. Then 'make' and you can see 'libxdrf.a'.
  * installed xdr file conversion supported by Gromacs
    * See [Gromacs Website](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library) if you want to get a recent version.
    * For convenience, you can get the v1.1.4 library, 'xdrfile-1.1.4.tar.gz', in **xdr** folder. Then install, for example:

```sh
$ ./configure --prefix=/home/(username)/xdr1.1.4
$ make install
```

After installing, add the following line to **.bashrc file** in your home directory.

```sh
export LD_LIBRARY_PATH=/home/(username)/xdr1.1.4/lib
```

  * Compiled with xtc normal or xtc wrapper version written by James Wes Barnett
    * In details, see ~~[Blog of James Wes Barnett](http://statthermo.blogspot.com)~~ broken currently. Instaed visit [GitHub Page](https://github.com/wesbarnett).
    * For convenience, the files are uploaded in **xtc** folder.
    * As for xtc normal version, How to compile is like: 

```sh
$ ifort -c xtc-interface.f90 -lxdrfile -L/home/hjung/xdr1.1.4/lib
```

    * To check, the command './a.out' with 'traj.xtc' file will print out Time, Step, Precision, No.Atoms for test. Then compile wrapper version of read_xtc_file after compiling:

```sh
ifort xtc-interface.o readxtc.f90 -lxdrfile -L/home/hjung/xdr1.1.4/lib
```

    * As for xtc wrapper version, compile like:

```sh
$ ifort -c xtc-interface-wrap.f90 -L/home/hjung/xdr1.1.4/lib -lxdrfile
```

      Test mode:
```sh
$ ifort xtc-interface-wrap.o readxtc-wrap.f90 -lxdrfile -L/home/hjung/xdr1.1.4/lib
```

      To test, do './a.out' and/or compare with previous result './a.out'.
    * Change of original code of wrapper version I did is that I modify module name *xtc* to *xtc_mod* in source files because intel fortran compiler prints errors that 'defined type, xtc is the same as module xtc.'. I think it comes from dulicate module names of two versions.
    *In my case, I choose wrapper version for xtc conversion. Either is fine.

## Install
It needs DensProf.f90, xtc-interface-wrap.o, and xtc_mod.mod (the last two file comes from compiling xtc code) to compile. 'xtc_mod.mod' file should be located in the same folder with DensProf.f90.

```sh
$ ifort -O2 DensProf.f90 xtc-interface-wrap.o -L/home/hjung/xdr1.1.4/lib -lxdrfile -o D-Prof.x 
```

## Usage for help
```sh
$ ./D-Prof.x
```


