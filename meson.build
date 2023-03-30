project('tiny_mc', 'c',
        default_options: [ 'c_std=c11', 'warning_level=2' ]
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m')

tiny_mc = executable('tiny_mc',
    sources: [ 'src/tiny_mc.c', 'src/wtime.c' ],
    include_directories: ['src'],
    dependencies: m_dep
)
