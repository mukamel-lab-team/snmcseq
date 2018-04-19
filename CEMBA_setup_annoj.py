#!/usr/bin/env python3
# Fangming Xie
# Caution: This script is indented by "2 tab" or "4 space" ? 


import subprocess as sp
from natsort import natsorted
# import os

from __init__ import *
from snmcseq_utils import cd
import snmcseq_utils
from snmcseq_utils import create_logger
from CEMBA_update_mysql import connect_sql


TEMPLATE_MC_PHP = (
"""<?php
$append_assembly = true;
$table = 'CEMBA_annoj.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("banjo.ucsd.edu", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_wgta.php';
?>
"""
)

TEMPLATE_BED_PHP = (
"""<?php
$append_assembly = true;
$table = 'CEMBA_annoj.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("banjo.ucsd.edu", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_masks.php';
?>
"""
)
TEMPLATE_DMR_PHP = TEMPLATE_BED_PHP
TEMPLATE_ATAC_PHP = TEMPLATE_BED_PHP

TEMPLATE_INDEX_MC = (
"""
    {{
      id   : '{}',
      name : '{}',
      type : 'MethTrack',
      path : 'DNA Methylation',
      data : './browser/fetchers/mc/{}.php',
      iconCls : 'salk_meth',
      height : 30,
      class : 'CG -CHH -CHG',
      showControls : true
    }},
"""
)

TEMPLATE_INDEX_DMR = (
"""
    {{
      id : '{}',
      name : '{}',
      type : 'ReadsTrack',
      path : 'DMR',
      data : './browser/fetchers/dmrs/{}.php',
      iconCls : 'salk_dmr',
      height : 5,
      scale : 100,
      color : {{read : '000099'}},
      showControls : 0,
      single : true,
    }},
"""
)

TEMPLATE_INDEX_ATAC = (
"""
      {{
        id : '{}',
        name : '{}',
        type : 'PairedEndTrack',
        path : 'ATAC-Seq',
        data : './browser/fetchers/atac/{}.php',
        iconCls : 'silk_bricks',
        height : 20,
        scale : 2,
        single : true,
      }},
"""
)


def gen_annoj_config(track_name, track_vis_name, template_string):
    """
    """
    return template_string.format(track_name, track_vis_name, track_name)

def gen_annoj_configs(track_names, track_vis_names, template_string):
    """
    """
    string = ''
    for track_name, track_vis_name in zip(track_names, track_vis_names):
        string += template_string.format(track_name, track_vis_name, track_name)
    return string

def gen_active_configs(lst):
    string = ''
    for item in lst:
        string +="\t'{}',\n".format(item)
    return string

def gen_php_template(track_name, template_string):
    """
    """
    return template_string.format(track_name, track_name, track_name)


TEMPLATE_INDEX_HTML=(
"""<html>
<head>
  <meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>
  <title>{}</title>


  <link type='text/css' rel='stylesheet' href='aj2/aj2_css_dev/ext-all.css' />
  <script type='text/javascript' src='aj2/js/ext-base-3.2.js'></script>
  <script type='text/javascript' src='aj2/js/ext-all-3.2.js'></script>

  <link type='text/css' rel='stylesheet' href='aj2/aj2_css_dev/viewport.css' />
  <link type='text/css' rel='stylesheet' href='aj2/aj2_css_dev/plugins.css' />
  <link type='text/css' rel='stylesheet' href='aj2/aj2_css_dev/salk.css' />
  <script type='text/javascript' src='aj2/js/excanvas.js'></script>
  <script type='text/javascript' src='aj2/js/aj-cndd1-dev.js'></script>

  <!-- Favicon -->
  <link rel="icon" href="/var/www/html/annoj/browser/aj.ico" type="image/x-icon">
  <link rel="shortcut icon" href="/var/www/html/annoj/browser/aj.ico" type="image/x-icon">


  <!-- Config -->
  <script type='text/javascript'>

  AnnoJ.config = {{

    info : {{
      title  : '{}',
      genome  : 'mm10',
      contact  : '',
      email : '',
      institution : ''
    }},
    tracks : [
      //Models
      {{
        id   : 'gene_model_mm10',
        name : 'Gene Models (mm10)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/genes_mm10.php',
        height : 100,
        showControls : true
      }},


      // DNA methylation
      {}
      // DMRs
      {}


    ],

    active : [// *** Gene models
      'gene_model_mm10',
      {}

    ],
    genome    : './browser/fetchers/mus_musculus.php',
    bookmarks : './browser/fetchers/mus_musculus.php',
    stylesheets : [
      {{
        id   : 'css1',
        name : 'Plugins CSS',
        href : 'css/plugins.css',
        active : true
      }},{{
        id   : 'css2',
        name : 'SALK CSS',
        href : 'css/salk.css',
        active : true
      }}
    ],
    location : {{
      assembly : '1',
      position : '132487763',
      bases    : 300,
      pixels   : 1
    }},
    admin : {{
      name  : 'Eran Mukamel',
      email : 'emukamel@ucsd.edu',
      notes : 'University of California, San Diego'
    }},
  }};



  </script>



</head>

<body>

  <noscript>
    <table id='noscript'><tr>
      <td><img src='hs/img/Anno-J.jpg' /></td>
      <td>
        <p>Anno-J cannot run because your browser is currently configured to block Javascript.</p>
        <p>To use the application please access your browser settings or preferences, turn Javascript support back on, and then refresh this page.</p>
        <p>Thankyou, and enjoy the application!<br /></p>
      </td>
    </tr></table>
  </noscript>

  <!-- Enable URL queries -->
  <script type='text/javascript'> var queryPost; </script>
  <script type='text/javascript' src='./browser/js/urlinit.js'></script>
  <!-- Google Analytics -->
  <script type="text/javascript">
  var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
  document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
  </script>

  <script type="text/javascript">
  var pageTracker = _gat._getTracker("UA-4150298-1");
  pageTracker._initData();
  pageTracker._trackPageview();
  </script>
  <script type="text/javascript">
  var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
  document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
  </script>
  <script type="text/javascript">
  try {{
    var pageTracker = _gat._getTracker("UA-4150298-1");
    pageTracker._trackPageview();
  }} catch(err) {{}}</script>
</body>

</html>
""" # .format(title, title, mc, dmr, active)
)


def gen_index_html(template_string, title, mc, dmr, active):
    return (template_string.format(title, title, mc, dmr, active))

def setup_annoj_main(ens):
    """Given an ensemble, setup annoj browser 
    """

    log = create_logger()
    time_stamp = time.time()

    logging.info("Setting up {} annoj browser...".format(ens))
    ### split allc tables of an ensemble to /scratch/tmp_allc_split
    # /scratch/tmp_allc_split needs to be removed
    output_folder = '/scratch/tmp_allc_split_{}'.format(time_stamp)
    try:
        os.makedirs(output_folder)
        logging.info('{} created!'.format(output_folder))
    except:
        raise ValueError

    logging.info("Splitting allc files to {}".format(output_folder))
    for allc_file in sorted(glob.glob(os.path.join(PATH_ENSEMBLES, ens, 'allc_merged/allc_merged_*.tsv'))):
        # split allc files
        logging.info('Processing {}...'.format(os.path.basename(allc_file)))
        df = snmcseq_utils.read_allc_CEMBA(allc_file, pindex=False, compression='infer')
        for chrom, df_chrom in df.groupby('chr'):
            output = '{}/{}_{}.tsv'.format(output_folder, 
                os.path.basename(allc_file)[:-len('.tsv')], chrom) 
            df_chrom.to_csv(output, sep="\t", na_rep='NA', index=False, header=False)

    # rename splitted allc for further processes (added clst removed merged)
    allc_split_files = sorted(glob.glob(os.path.join(output_folder, 'allc_merged_*.tsv')))
    for file in allc_split_files:
        src = file
        file = os.path.basename(file)
        file_l = file.split('_')
        file_l[-3] = 'clst' + file_l[-3]
        file = '_'.join(file_l)
        dst = os.path.join(output_folder, file)
        os.rename(src, dst)
    allc_split_files = sorted(glob.glob(os.path.join(output_folder, 'allc_merged_*.tsv')))


    ### allc to mc
    # mc files are in /scratch -> /scratch/tmp_mc_split_$time
    logging.info("Allc to mc files...")
    allc_split_files = sorted(glob.glob(os.path.join(output_folder, 'allc_merged_*.tsv')))
    mc_split_files = ['mc_' + os.path.basename(file)[len('allc_'):] 
                      for file in allc_split_files]

    with cd(output_folder):
        cmd = '/cndd/fangming/lab_scripts/allc2mc.pl ./allc_merged_*.tsv'
        sp.call(cmd, shell=True)

    # move all relevent mc files to a tmp folder
    output_folder_mc = '/scratch/tmp_mc_{}'.format(time_stamp) 
    try:
        os.makedirs(output_folder_mc)
        logging.info('{} created!'.format(output_folder_mc))
    except:
        raise ValueError
        
    for mc_split_file in mc_split_files:
        src = os.path.join('/scratch', mc_split_file)
        dst = os.path.join(output_folder_mc, mc_split_file)
        os.rename(src, dst)

    ### mc to mysql 
    logging.info('Mc files to mySQL...')
    with cd(output_folder_mc):
        cmd = '/cndd/fangming/lab_scripts/load_mysql_mc CEMBA_annoj > /cndd/fangming/CEMBA/snmcseq_dev/logs/annoj_1.log'
        sp.call(cmd, shell=True)

    # remove mc tmp files
    cmd = 'rm -r {}'.format(output_folder_mc)
    sp.call(cmd, shell=True)

    ### dmr to bed
    num_dms = NUM_DMS
    logging.info("DMR to bed files...")
    cluster_type = 'cluster_mCHmCG_lv_npc50_k30'
    annotation_type = 'annotation' + cluster_type[len('cluster'):]
    engine = connect_sql('CEMBA')
    output_folder_dmr = '/scratch/tmp_dmr_{}'.format(time_stamp)

    # cluster annotation
    sql = """SELECT {}, {} FROM {}""".format(cluster_type, annotation_type, ens)
    df_info = pd.read_sql(sql, engine)
    df_info = df_info.sort_values(cluster_type).drop_duplicates().set_index(cluster_type)

    ens_path = os.path.join(PATH_ENSEMBLES, ens)
    n_clusters = len(glob.glob(os.path.join(ens_path, 'allc_merged/allc_merged_mCG_{}_*_{}.tsv'.format(cluster_type, ens))))
    # dmr results
    input_f = os.path.join(ens_path, 'dmr/dmr_allc_merged_mCG_{}_rms_results_collapsed.tsv'.format(cluster_type))
    df = pd.read_table(input_f, index_col=['#chr', 'start', 'end'], dtype={'#chr': object})
    # df_hypo
    df_hypo = df.loc[(df['number_of_dms']>=num_dms) & (~df['hypomethylated_samples'].isnull()), 'hypomethylated_samples'].apply(
            lambda x: x.split(','))
    try:
        os.makedirs(output_folder_dmr)
        logging.info('Created directory: {}'.format(output_folder_dmr))
    except:
        raise ValueError

    for i in range(n_clusters):
        df_hypo_cluster = df_hypo.loc[
            df_hypo.apply(lambda x: ('merged_mCG_{}_{}_{}'.format(cluster_type, i+1, ens) in x))]
        output = os.path.join(output_folder_dmr, 'dmr_hypo_allc_merged_mCG_cluster_mCHmCG_lv_npc50_k30_clst{}.bed'.format(i+1))
        df_hypo_cluster.to_csv(output, sep='\t', header=False, index=True, na_rep='NA')
        logging.info("Saved to {}".format(output))
        logging.info(df_hypo_cluster.shape)

    ### rename bed
    with cd(output_folder_dmr):
        for file in glob.glob('*'):
            src = file 
            dst = '{}_{}.bed'.format(src[:-len('.bed')], ens).replace('allc_merged_', '') 
            os.rename(src, dst)

    ### bed to mysql
    logging.info("Uploading bed files (DMR) to mySQL databases...")
    bed_files = sorted(glob.glob(os.path.join(output_folder_dmr, '*')))
    with cd(output_folder_dmr):
        for bed_file in bed_files:
            # fn, db, table_prefix
            cmd = '/cndd/fangming/lab_scripts/bed_load_mysql.sh {} {} {}'.format(os.path.basename(bed_file), 
                                                                                 'CEMBA_annoj', 
                                                                                 os.path.basename(bed_file)[:-len('.bed')])
            sp.call(cmd, shell=True)

    # remove mc tmp files
    cmd = 'rm -r {}'.format(output_folder_dmr)
    sp.call(cmd, shell=True)


    ### create php file
    path = '/var/www/html/annoj_private/CEMBA'

    # php file for mc_split_file
    output_folder_php = '/cndd/fangming/CEMBA/snmcseq_dev/browsers/php_mc_{}'.format(ens)
    if not os.path.isdir(output_folder_php):
        os.makedirs(output_folder_php)
        logging.info('Created directory: {}'.format(output_folder_php))
        
    track_names_mc = np.unique(['_'.join(os.path.basename(file).split('_')[:-1]) for file in mc_split_files])
    for track_name in track_names_mc:
        with open(os.path.join(output_folder_php, '{}.php'.format(track_name)), 'w') as f_php:
            f_php.write(gen_php_template(track_name, TEMPLATE_MC_PHP))

    # php file for dmr files
    output_folder_php = '/cndd/fangming/CEMBA/snmcseq_dev/browsers/php_dmrs_{}'.format(ens)
    if not os.path.isdir(output_folder_php):
        os.makedirs(output_folder_php)
        logging.info('Created directory: {}'.format(output_folder_php))
    track_names_dmr = [os.path.basename(file)[:-len('.bed')] for file in bed_files]
    for track_name in track_names_dmr:
        with open(os.path.join(output_folder_php, '{}.php'.format(track_name)), 'w') as f_php:
            f_php.write(gen_php_template(track_name, TEMPLATE_DMR_PHP))

    ### index_file
    track_vis_names_mc = [file.replace('merged_', '').replace('mc_', '') for file in track_names_mc]
    track_vis_names_dmr = track_names_dmr 

    title = 'CEMBA browser ({})'.format(ens)
    mc_configs = gen_annoj_configs(track_names_mc, track_vis_names_mc, TEMPLATE_INDEX_MC) 
    dmr_configs = gen_annoj_configs(track_names_dmr, track_vis_names_dmr, TEMPLATE_INDEX_DMR) 

    # active config, active order (natsorted mc track and dmr tracks)
    active_orders = np.vstack((np.asarray(natsorted(track_names_mc)), 
                               np.asarray(natsorted(track_names_dmr)))).flatten('F')
    active_configs = gen_active_configs(active_orders) 

    # gen_index
    index_html = gen_index_html(TEMPLATE_INDEX_HTML, title, mc_configs, dmr_configs, active_configs)

    # output 
    output_html = '/cndd/fangming/CEMBA/snmcseq_dev/browsers/index_{}.html'.format(ens)
    with open(output_html, 'w') as f:
        f.write(index_html)

    logging.info("Setting up {} annoj browser is done!".format(ens))

if __name__ == '__main__':

    enss = ['Ens8', 'Ens9', 'Ens10']
    for ens in enss:
        setup_annoj_main(ens)
