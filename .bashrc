# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

login_test=$(echo ${HOSTNAME} | grep login | wc -l);
if [[ $login_test == 1 ]];then 
    qinteract;
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
# Adding local version of perlbrew:

#perl_parent=/home/jas297/linux/;
perl_parent=/mnt/clustertmp/common/rja20_dev/;
# SAMBA support

PERLBREW_HOME=${perl_parent}/perl5;
PERLBREW_ROOT=${perl_parent}/perl5;

export PERLBREW_HOME;
export PERLBREW_ROOT;

source "${perl_parent}/perl5/etc/bashrc";

PERL_MB_OPT="--install_base ${perl_parent}/perl5"; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=${perl_parent}/perl5"; export PERL_MM_OPT;
PERL5LIB="${perl_parent}/perl5/lib/perl5"; export PERL5LIB;
PATH="${perl_parent}/perl5/bin:$PATH"; export PATH;
PERL_LOCAL_LIB_ROOT="${perl_parent}/perl5:$PERL_LOCAL_LIB_ROOT"; export PERL_LOCAL_LIB_ROOT;


#eval "$(perl -I${PERL5LIB} -Mlocal::lib)"


# Updating env variable for MATLAB2015b runtime:
codepath=/mnt/clustertmp/common/rja20_dev/;
#export LD_LIBRARY_PATH=${codepath}/MATLAB2015b_runtime/v90/runtime/glnxa64:${codepath}/MATLAB2015b_runtime/v90/bin/glnxa64:${codepath}/MATLAB2015b_runtime/v90/sys/os/glnxa64:${LD_LIBRARY_PATH};
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${codepath}/MATLAB2015b_runtime/v90/runtime/glnxa64:${codepath}/MATLAB2015b_runtime/v90/sys/os/glnxa64;
# Adding local version of perlbrew:
#source ./perl5/etc/bashrc

# Adding local version of perlbrew:
source "${perl_parent}/perl5/etc/bashrc";


export SAMBA_PATH=/mnt/clustertmp/common/rja20_dev/SAMBA/;
source "${SAMBA_PATH}/bashrc_for_SAMBA";
