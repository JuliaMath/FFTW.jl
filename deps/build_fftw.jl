using Libdl
using BinDeps
using BinDeps: builddir, depsdir, libdir

# Binaries is not a recognized provider on Linux >:/
modified_defaults = false
if !in(BinDeps.Binaries, BinDeps.defaults)
    pushfirst!(BinDeps.defaults, BinDeps.Binaries)
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
    nover = replace(lib, major => "")
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
                              "43a049496d47dd1919283a981d21078942d04b067edd75b61cc6aef538098b53"),
    "i686-pc-linux-gnu"   => ("$URL-linux-i686.tar.gz",
                              "0b6ec1440ea76f72ed6e5d105a4d634e195b484e847c521ac0812092069f0c74"),
    "x86_64-apple-darwin" => ("$URL-osx-x86_64.tar.gz",
                              "193f20bba844fd2524f70e941aa5203b153bdcd321df2f20fc3a9d15b4aae50a"),
    "x86_64-w64-mingw32"  => ("$URL-win-x86_64.zip",
                              "4e9bd1022c5980869b06bb9d10767fd2578dfe47f9323e6f572d0918da8c8a07"),
    "i686-w64-mingw32"    => ("$URL-win-i686.zip",
                              "60a65fd689495629a1558858f3a52011f9ce28df8aa756f5b315b17dd8dd8843"),
)

const machine = Sys.isapple() ? "x86_64-apple-darwin" : Sys.MACHINE

if haskey(downloads, machine)
    url, sha = downloads[machine]
    provides(Binaries, URI(url), [libfftw, libfftwf], SHA=sha, os=BinDeps.OSNAME,
             unpacked_dir=joinpath("usr", "lib"), installed_libpath=libdir(libfftw))
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
    popfirst!(BinDeps.defaults)
end
