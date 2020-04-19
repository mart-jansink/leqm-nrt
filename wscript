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

def configure(conf):
    conf.load('compiler_cxx')
    conf.env.append_value('CXXFLAGS', ['-Wall', '-Wextra', '-D_FILE_OFFSET_BITS=64', '-D__STDC_FORMAT_MACROS'])
    conf.env.append_value('LINKFLAGS', ['-pthread'])
    conf.check_cfg(package='sndfile', args='--cflags --libs', uselib_store='SNDFILE')

def build(bld):
    bld(source='libleqm_nrt.pc.in',
        version=VERSION,
        includedir='%s/include/libdcp%s' % (bld.env.PREFIX, bld.env.API_VERSION),
        libs="-L${libdir} -lleqm_nrt",
        install_path='${LIBDIR}/pkgconfig')

    bld.recurse('src')
