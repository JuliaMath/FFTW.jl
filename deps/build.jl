using BinDeps
using BinDeps: builddir, usrdir

# Binaries is not a recognized provider on Linux >:/
modified_defaults = false
if !in(BinDeps.Binaries, BinDeps.defaults)
    unshift!(BinDeps.defaults, BinDeps.Binaries)
    modified_defaults = true
end

BinDeps.@setup

const FFTW_VER = v"3.3.6-pl2"

if Sys.iswindows()
    const libfftw_name = "libfftw3"
    const libfftwf_name = "libfftw3f"
else
    const libfftw_name = "libfftw3_threads"
    const libfftwf_name = "libfftw3f_threads"
end

# Why can't everyone just agree on what to call this library...
function makealiases(lib)
    major = string(FFTW_VER.major)
    nover = replace(lib, major, "")
    return String[
        nover,
        join([lib, Libdl.dlext, major], "."),
        join([nover, Libdl.dlext, major], "."),
        lib * "-" * major,
        nover * "-" * major,
    ]
end

libfftw = library_dependency(libfftw_name, aliases=makealiases(libfftw_name))
libfftwf = library_dependency(libfftwf_name, aliases=makealiases(libfftwf_name))

const URL = "https://github.com/ararslan/fftw-builder/releases/download/v$FFTW_VER/libfftw-$FFTW_VER"

# Mapping of Sys.MACHINE to (url, sha) for precompiled binaries from fftw-builder
const downloads = Dict(
    "x86_64-pc-linux-gnu" => ("$URL-linux-x86_64.tar.gz",
                              "b0576a403691a54bb3c234cccb8f0b97bb8f51ca1835c43c69e9c90358eb2647"),
    "i686-pc-linux-gnu"   => ("$URL-linux-i686.tar.gz",
                              "d7ce5f18719eb93691ba5c3f7c04b97dde54e4569581a063be096ea52f8126e5"),
    "x86_64-apple-darwin" => ("$URL-osx-x86_64.tar.gz",
                              "6f354ed992e9c8f33b565d6bc32ea0ff08c8caf8e4a00402c0c1ab56106c706e"),
    "x86_64-w64-mingw32"  => ("$URL-win-x86_64.zip",
                              "915a0290d495d8d1c040bec8743ccfe35246fe85bdba79ac11b66f2f6a413ab3"),
    "i686-w64-mingw32"    => ("$URL-win-i686.zip",
                              "84ecb82afb9f210a711391aad11b76af7b8f284d7e56f5d5d409e1c47cc547c2"),
)

const machine = Sys.isapple() ? "x86_64-apple-darwin" : Sys.MACHINE

if haskey(downloads, machine)
    url, sha = downloads[machine]
    isdir(usrdir(libfftw)) || mkpath(usrdir(libfftw))
    provides(Binaries, URI(url), [libfftw, libfftwf], SHA=sha, os=BinDeps.OSNAME,
             unpacked_dir="fftw-$FFTW_VER")
    scratch = false
elseif Sys.KERNEL === :FreeBSD
    provides(BSDPkg, "fftw3", [libfftw, libfftwf], os=:FreeBSD)
    scratch = false
else
    info("No precompiled binaries found for your system. Building from scratch...")
    scratch = true
end

general_config = ["--prefix=" * abspath(builddir(libfftw)),
                  "--libdir=" * abspath(libdir(libfftw)),
                  "--bindir=" * abspath(bindir(libfftw))]

fftw_config = ["--enable-shared", "--disable-fortran", "--disable-mpi", "--enable-threads"]
fftw_enable_single = "--enable-single"

if Sys.ARCH === :ppc
    append!(fftw_config, ["--enable-altivec", "--enable-fma"])
elseif Sys.ARCH === :x86_64
    append!(fftw_config, ["--enable-sse2", "--enable-fma"])
end

if Sys.iswindows()
    append!(fftw_config, ["--with-our-malloc", "--with-combined-threads"])
    Sys.ARCH === :x86_64 || push!(fftw_config, "--with-incoming-stack-boundary=2")
end

# Make it harder to build from scratch
if scratch
    provides(Sources, URI("http://www.fftw.org/fftw-$FFTW_VER.tar.gz"), [libfftw, libfftwf])

    provides(BuildProcess, (@build_steps begin
        GetSources(libfftw)
        CreateDirectory(builddir(libfftw))
        @build_steps begin
            ChangeDirectory(builddir(libfftw))
            FileRule(joinpath(libdir(libfftw), libfftw_name * "." * Libdl.dlext), @build_steps begin
                CreateDirectory(libdir(libfftw))
                `$(joinpath(srcdir(libfftw), "fftw-$FFTW_VER", "configure")) $general_config $fftw_config`
                `$MAKE_CMD`
                `$MAKE_CMD install`
            end)
            FileRule(joinpath(libdir(libfftw), libfftwf_name * "." * Libdl.dlext), @build_steps begin
                `$(joinpath(srcdir(libfftw), "fftw-$FFTW_VER", "configure")) $general_config $fftw_config $fftw_enable_single`
                `$MAKE_CMD`
                `$MAKE_CMD install`
            end)
        end
    end), [libfftw, libfftwf])
end

if Sys.iswindows()
    BinDeps.@install Dict([:libfftw3 => :libfftw, :libfftw3f => :libfftwf])
else
    BinDeps.@install Dict([:libfftw3_threads => :libfftw, :libfftw3f_threads => :libfftwf])
end

if modified_defaults
    shift!(BinDeps.defaults)
end
