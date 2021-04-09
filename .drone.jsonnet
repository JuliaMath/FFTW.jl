local Pipeline(os, arch, version) = {
    kind: "pipeline",
    name: os+" - "+arch+" - Julia "+version,
    platform: {
	os: os,
	arch: arch
    },
    steps: [
	{
	    name: "build",
	    image: "julia:"+version,
	    commands: [
		"julia --project=. --check-bounds=yes --color=yes -e 'using InteractiveUtils; versioninfo(verbose=true); using Pkg; Pkg.test(coverage=true)'"
	    ]
	}
    ],
    trigger: {
        branch: ["master"]
    }
};

[
    # Commenting this out because we don't have an official armv7l build yet
    #Pipeline("linux", "arm",   "1.6"),
    Pipeline("linux", "arm64", "1.6")
]
