from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="ntdmc_trachoma.libtrachoma",
            sources=[
                "src/libtrachoma/periods.c",
                "src/libtrachoma/shift.c",
                "src/libtrachoma/trachoma.c"
            ],
        ),
    ]
)
