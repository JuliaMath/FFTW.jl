using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libfftw3"], :libfftw3),
    LibraryProduct(prefix, ["libfftw3f"], :libfftw3f),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaBinaryWrappers/FFTW_jll.jl/releases/download/FFTW-v3.3.9+2"

# Listing of files generated by BinaryBuilder:
download_info = Dict(
    Linux(:aarch64, libc=:glibc) => ("$bin_prefix/FFTW.v3.3.9.aarch64-linux-gnu.tar.gz", "ca6ed03d45a1b4f054c32d884ab3a20f55bd20c166372a35ceefe91f2f10c94d"),
    Linux(:aarch64, libc=:musl) => ("$bin_prefix/FFTW.v3.3.9.aarch64-linux-musl.tar.gz", "7d6400a4df54b671f301b80b9e79bee7875cfd4a78e519fa6dbfdd38b4ccfc28"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf) => ("$bin_prefix/FFTW.v3.3.9.arm-linux-gnueabihf.tar.gz", "3c3f630c265e1806422bfe0bb5beebfe401651bf5c08ff2dd160cc859d3faf0d"),
    Linux(:armv7l, libc=:musl, call_abi=:eabihf) => ("$bin_prefix/FFTW.v3.3.9.arm-linux-musleabihf.tar.gz", "2b286cd2afcd68caf1dd64b1c9c333c124099c054c4903e4e3f68ad3d176f50a"),
    Linux(:i686, libc=:glibc) => ("$bin_prefix/FFTW.v3.3.9.i686-linux-gnu.tar.gz", "aa5e90cfc15a532b7cd6a1fa9ea79a6015c1d4044cf2c2f3d927523b2e8e5f56"),
    Linux(:i686, libc=:musl) => ("$bin_prefix/FFTW.v3.3.9.i686-linux-musl.tar.gz", "0d0ef98bbf1ff1955414dab22ef11a3163994baf3ffcfe9e8e9ac14dd141fda8"),
    Windows(:i686) => ("$bin_prefix/FFTW.v3.3.9.i686-w64-mingw32.tar.gz", "5b70660701b6466e0cbd3dd7f4fb741e2098110114cac8386f1c8dcd7b9c1b7f"),
    Linux(:powerpc64le, libc=:glibc) => ("$bin_prefix/FFTW.v3.3.9.powerpc64le-linux-gnu.tar.gz", "6a3e1d1cc34b24b329a91d0907ae4b50e8088aa83648ed251967c4a72f7b5390"),
    MacOS(:x86_64) => ("$bin_prefix/FFTW.v3.3.9.x86_64-apple-darwin14.tar.gz", "a1779b7c4ded4840544febb093bfb03da4b6ade37ac996979a48fd11d8c494a6"),
    Linux(:x86_64, libc=:glibc) => ("$bin_prefix/FFTW.v3.3.9.x86_64-linux-gnu.tar.gz", "5a7479ec03cb729a3cd51e9fc222d0bd478f54438d27d649a1091438649be4ca"),
    Linux(:x86_64, libc=:musl) => ("$bin_prefix/FFTW.v3.3.9.x86_64-linux-musl.tar.gz", "8c7acde8c76726a14db64ffc41357fa8d1950ba992ef2aa34ac08c462f7708b7"),
    FreeBSD(:x86_64) => ("$bin_prefix/FFTW.v3.3.9.x86_64-unknown-freebsd11.1.tar.gz", "26812b0ecda157a674856f1f671a5faef4bc9133a904fc1dca36bb32bcd6fad5"),
    Windows(:x86_64) => ("$bin_prefix/FFTW.v3.3.9.x86_64-w64-mingw32.tar.gz", "fdbb0002ff49458d9faec14488d01ff904ec003a0c937649972ba92dc72e9f04"),
)

# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)
dl_info = choose_download(download_info, platform_key_abi())
if dl_info === nothing && unsatisfied
    # If we don't have a compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
end

# If we have a download, and we are unsatisfied (or the version we're
# trying to install is not itself installed) then load it up!
if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
    # Download and install binaries
    install(dl_info...; prefix=prefix, force=true, verbose=verbose)
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)