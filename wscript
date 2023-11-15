import shlex
import subprocess

this_version = subprocess.Popen(shlex.split('git tag -l --points-at HEAD'), stdout=subprocess.PIPE).communicate()[0]
last_version = subprocess.Popen(shlex.split('git describe --tags --abbrev=0'), stdout=subprocess.PIPE).communicate()[0]

if this_version == '':
    VERSION = '%sdevel' % last_version[1:].strip()
else:
    VERSION = this_version[1:].strip()

def options(opt):
    opt.load('compiler_cxx')
    opt.add_option('--static', action='store_true', default=False, help='build libleqm_nrt statically')
    opt.add_option('--without-libsndfile', action='store_true', default=False, help='do not build code that requires libsndfile (file-based interface and CLI tool)')

def configure(conf):
    conf.load('compiler_cxx')
    conf.env.WITH_LIBSNDFILE = not conf.options.without_libsndfile
    conf.env.STATIC = conf.options.static
    if conf.env.WITH_LIBSNDFILE:
        conf.env.append_value('CXXFLAGS', ['-DLEQM_NRT_WITH_LIBSNDFILE'])
        conf.check_cfg(package='sndfile', args='--cflags --libs', uselib_store='SNDFILE')
    conf.env.append_value('CXXFLAGS', ['-Wall', '-Wextra', '-D_FILE_OFFSET_BITS=64', '-D__STDC_FORMAT_MACROS', '-std=c++11', '-O2'])
    conf.env.append_value('LINKFLAGS', ['-pthread'])

def build(bld):
    bld(source='leqm_nrt.pc.in',
        version=VERSION,
        includedir='%s/include' % bld.env.PREFIX,
        libs="-L${libdir} -lleqm_nrt",
        install_path='${LIBDIR}/pkgconfig')

    bld.recurse('src')
