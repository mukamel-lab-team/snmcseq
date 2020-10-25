"""AnnoJ browser templates
Fangming Xie

Available templates:

php:
- TEMPLATE_MC_PHP
- TEMPLATE_MC_SINGLE_PHP
- TEMPLATE_BED_PHP
- TEMPLATE_DMR_PHP
- TEMPLATE_ATAC_PHP
- TEMPLATE_RNA_PHP

html:
- TEMPLATE_INDEX_MC
- TEMPLATE_INDEX_MC_SINGLE
- TEMPLATE_INDEX_DMR
- TEMPLATE_INDEX_ATAC
- TEMPLATE_INDEX_RNA

- TEMPLATE_INDEX_HTML (mouse)
"""

import collections
import os

# Php file templates
TEMPLATE_MC_PHP = (
"""<?php
$append_assembly = true;
$table = '{}.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("{}", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_wgta.php';
?>
"""
)

TEMPLATE_MC_SINGLE_PHP = (
"""<?php
$append_assembly = true;
$table = '{}.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("{}", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_wgta_single.php';
?>
"""
)

TEMPLATE_BED_PHP = (
"""<?php
$append_assembly = true;
$table = '{}.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("{}", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_masks.php';
?>
"""
)

TEMPLATE_BEDGRAPH_PHP = (
"""<?php
$append_assembly = true;
$table = '{}.{}_';
$title = '{}';
$info = '{}';
$link = mysql_connect("{}", "cndd_annoj", "jonna_ddnc") or die("failed");
require_once '../includes/common_bedgraph.php';
?>
"""
)

TEMPLATE_DMR_PHP = TEMPLATE_BED_PHP
TEMPLATE_ATAC_PHP = TEMPLATE_BED_PHP
# TEMPLATE_RNA_PHP = TEMPLATE_BED_PHP
TEMPLATE_RNA_PHP = TEMPLATE_BEDGRAPH_PHP
TEMPLATE_DATABASE = 'CEMBA_annoj'

def gen_php_template(track_name, template_string, database="localhost", template_database=TEMPLATE_DATABASE):
    """
    """
    return template_string.format(template_database, track_name, track_name, track_name, database)


# Index block templates 
TEMPLATE_INDEX_MCG = (
"""
    {{
      id   : '{0}',
      name : '{1}',
      type : 'MethTrack',
      path : 'DNA Methylation',
      data : './browser/fetchers/mc/{2}.php',
      iconCls : 'salk_meth',
      height : 30,
      class : 'CG -CH',
      single: true,
      color : {{CG: '{3}', CH: '{3}'}},
      showControls: true,
    }},
"""
)
TEMPLATE_INDEX_MC = TEMPLATE_INDEX_MCG
TEMPLATE_INDEX_MCH = (
"""
    {{
      id   : '{0}',
      name : '{1}',
      type : 'MethTrack',
      path : 'DNA Methylation',
      data : './browser/fetchers/mc/{2}.php',
      iconCls : 'salk_meth',
      height : 30,
      class : 'CH -CG',
      single: true,
      color : {{CH: '{3}', CG: '{3}'}},
      showControls: true,
    }},
"""
)

# TEMPLATE_INDEX_MC_SINGLE = (
# """
#     {{
#       id   : '{}',
#       name : '{}',
#       type : 'MethTrack',
#       path : 'DNA Methylation',
#       data : './browser/fetchers/mc_single/{}.php',
#       iconCls : 'salk_meth',
#       height : 30,
#       class : 'CGdmr CG COV', 
#       single : true,
#       showControls : true
#     }},
# """
# )

TEMPLATE_INDEX_DMR = (
"""
    {{
      id : '{}',
      name : '{}',
      type : 'ReadsTrack',
      path : 'DMR',
      data : './browser/fetchers/dmr/{}.php',
      iconCls : 'salk_dmr',
      height : 5,
      scale : 100,
      color : {{read : '{}'}},
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
        color : {{count: '{}'}},
        scale : {},
        single : true,
      }},
"""
)

TEMPLATE_INDEX_RNA = (
"""
      {{
        id : '{}',
        name : '{}',
        type : 'PairedEndTrack',
        path : 'RNA-Seq',
        data : './browser/fetchers/rna/{}.php',
        iconCls : 'silk_bricks',
        height : 30,
        color : {{count: '{}'}},
        scale : {},
        single : true,
      }},
"""
)

TEMPLATE_INDEX_MODEL_HUMAN = (
"""
      {
        id   : 'gene_model_hg19',
        name : 'Gene Models (hg19)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/genes_hg19.php',
        height : 100,
        showControls : true
      },
      {
        id   : 'gene_model_hg19',
        id   : 'models_hg19_gencode19_gene',
        name : 'Gene Models (hg19 Gencode)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/models_hg19_gencode19_gene.php',
        height : 100,
        showControls : true
      },
      {
        id   : 'models_hg19_gencode_tpx',
        name : 'Transcript Models (hg19 Gencode)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/models_hg19_gencode_tpx.php',
        height : 150,
        showControls : true
      },
      {
        id   : 'CGi',
        name : 'CpG islands',
        type : 'ReadsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/hg19_CGi.php',
        height : 20,
        scale : 100,
        showControls : true
      },
      {
        id   : 'Repeats',
        name : 'Repeats',
        type : 'ReadsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/hg19_repeats.php',
        height : 20,
        scale : 100,
        showControls : true
      },
      {
        id   : 'gene_model_hg19',
        name : 'Gene Models (hg19)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/genes_hg19.php',
        height : 100,
        showControls : true
      },
"""
)

TEMPLATE_INDEX_MODEL_MOUSE = (
"""
      {
        id   : 'gene_model_mm10',
        name : 'Gene Models (mm10)',
        type : 'ModelsTrack',
        path : 'Annotation models',
        data : './browser/fetchers/models/genes_mm10.php',
        height : 100,
        showControls : true
      },
"""
)

def gen_annoj_config(track_name, track_vis_name, track_color, template_string, use_scale=False, scale=1):
    """
    - scale: ingored if use_scale=False
    """
    if use_scale:
        return template_string.format(track_name, track_vis_name, track_name, track_color, scale)
    else:
        return template_string.format(track_name, track_vis_name, track_name, track_color)

def gen_annoj_configs(track_names, track_vis_names, track_colors, template_string, use_scale=False, scales='default', comment='Tracks'):
    """
    Args:
     - scales: a num or a list of numbers 
    """
    if isinstance(scales, str) and scales == 'default':
    	scales = [1]*len(track_names)

    string = '      //{}\n'.format(comment)
    for track_name, track_vis_name, track_color, scale in zip(track_names, track_vis_names, track_colors, scales):
        string += gen_annoj_config(track_name, track_vis_name, track_color, template_string, use_scale=use_scale, scale=scale)
    return string

def gen_active_configs(lst):
    """
    """
    string = ''
    for item in lst:
        string +="\t'{}',\n".format(item)
    return string

TEMPLATE_INDEX_HTML=(
"""<html>
<head>
  <meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>
  <title>{0}</title>

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
      title  : '{0}',
      genome  : '',
      contact  : '',
      email : '',
      institution : ''
    }},
    tracks : [
      //Models
      {3}


      //tracks
      {1}

    ],

    active : [// *** Gene models
      'gene_model_mm10',
      {2}

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
""" # .format(title, body, active, model)
)

TEMPLATE_INDEX_HTML_MOUSE = TEMPLATE_INDEX_HTML

TEMPLATE_INDEX_HTML_HUMAN=(
"""<html>
<head>
  <meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>
  <title>{0}</title>


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
      title  : '{0}',
      genome  : '',
      contact  : '',
      email : '',
      institution : ''
    }},
    tracks : [
      //Models
      {3}

      //tracks
      {1}

    ],

    active : [// *** Gene models
      'models_hg19_gencode_tpx',
      {2}

    ],
    genome    : './browser/fetchers/homo_sapiens_hg19.php',
    bookmarks : './browser/fetchers/homo_sapiens_hg19.php',
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
      assembly : '12',
      position : '113899838',
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
""" # .format(title, body, active, model)
)





def gen_index_html(template_string, title, body, active, template_model):
    """
        Args:
            - body: a string generated from gen_annoj_configs 
    """
    return (template_string.format(title, body, active, template_model))

# end bricks


TrackInfoT = collections.namedtuple('TrackInfoT', 
    [
    'track_comment',
    'php_template',
    'index_template',
    'track_names', 
    'track_vis_names', 
    'track_colors',
    'scales', 
    'use_scale',
    ])

def TrackInfo(php_template, index_template, 
    track_names, track_vis_names, track_colors, 
    scales='default', use_scale=False, track_comment='Track'):
    """Make TrackInfoT namedtuple
    Args:
        - track_names 
        - track_vis_names 
        - scales: 'default' or a list of numbers
        - use_scale
    """
    assert len(track_names) == len(track_vis_names)

    if use_scale and not isinstance(scales, str):
        assert len(scales) == len(track_names)
    else:
        scales = [1]*len(track_names)

    track_info = TrackInfoT(**{
        'php_template': php_template,
        'index_template': index_template,
        'track_names': track_names, 
        'track_vis_names': track_vis_names, 
        'track_colors': track_colors, 
        'scales': scales, 
        'use_scale': use_scale, 
        'track_comment': track_comment,
        })
    return track_info


def setup_annoj_worker(title, track_info_list, output_php_dir, output_html, 
    active_track_list=[], database='localhost', 
    template_html=TEMPLATE_INDEX_HTML_MOUSE,
    template_database=TEMPLATE_DATABASE, 
    template_model=TEMPLATE_INDEX_MODEL_MOUSE,
    ):
    """
    Args:
        - track_info_dict: a dictionary contains many "Track_info" namedtuple
        - active_tracks: a list 

    It generates
    """
    # check if every value is a TrackInfoT 
    for track_info in track_info_list:
        if not isinstance(track_info, TrackInfoT):
            raise ValueError("This is not a TrackInfo namedtuple.")

    # gen phps
    if not os.path.isdir(output_php_dir):
        os.makedirs(output_php_dir)

    for track_info in track_info_list:
        # one track type
        if not os.path.isdir(os.path.join(output_php_dir, 'php_'+track_info.track_comment)):
            os.makedirs(os.path.join(output_php_dir, 'php_'+track_info.track_comment))

        for track_name in track_info.track_names:
            # one track
            output_php_file = os.path.join(
                output_php_dir, 
                'php_'+track_info.track_comment, 
                track_name+'.php',
                )
            with open(output_php_file, 'w') as fphp:
                fphp.write(gen_php_template(
                    track_name, 
                    track_info.php_template, 
                    database=database, 
                    template_database=template_database,
                    ))

    # gen htmls
    body = ''
    for track_info in track_info_list: 
        body += gen_annoj_configs(
            track_info.track_names, 
            track_info.track_vis_names, 
            track_info.track_colors, 
            track_info.index_template, 
            use_scale=track_info.use_scale,
            scales=track_info.scales,
            comment=track_info.track_comment,
            )

    if not list(active_track_list):
        active_tracks = ''
    else:
        active_tracks = gen_active_configs(active_track_list)

    with open(output_html, 'w') as fhtml:
        fhtml.write(gen_index_html(template_html, title, body, active_tracks, template_model))

    return 
