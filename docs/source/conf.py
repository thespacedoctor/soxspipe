# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import time
import os
import sys
from datetime import datetime, date
# Check if we are on Read the Docs
on_rtd = os.environ.get('READTHEDOCS') == 'True'

# WHERE DOES THIS conf.py FILE LIVE?
moduleDirectory = os.path.dirname(os.path.realpath(__file__))
# GET PACKAGE __version__ INTO locals()
exec(open(moduleDirectory + "/../../soxspipe/__version__.py").read())
sys.path.insert(0, os.path.abspath('../../soxspipe/soxspipe'))

# General information about the project.
now = datetime.now()
now = now.strftime("%Y")
project = u'soxspipe'
copyright = u'%(now)s, David R. Young & Marco Lanodi' % locals()
version = "v" + str(__version__)
release = version
today_fmt = '%Y'

# -- General configuration ---------------------------------------------------
templates_path = ['_templates']
exclude_patterns = ['**xxx**']
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

extensions = [
    # Sphinx's own extensions
    "sphinx.ext.autodoc",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    'sphinxcontrib.inkscapeconverter',
    'sphinx.ext.autosectionlabel',
    # External stuff
    "sphinxext.opengraph",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_inline_tabs",
    'autodoc2',
    'sphinxcontrib.mermaid',
    'sphinx_togglebutton',
    'sphinx.ext.coverage',
    'sphinx.ext.linkcode',
    'sphinx_search.extension',
    'sphinx_tippy',
    "sphinx_remove_toctrees",
    'sphinxcontrib.bibtex'
]
myst_enable_extensions = [
    "tasklist",
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "attrs_block"

]
suppress_warnings = ["myst.strikethrough"]
autosectionlabel_prefix_document = True

numfig = True

numfig_format = {
    'code-block': 'Listing %s',
    'figure': 'Fig. %s',
    'table': 'Table %s',
    'section': 'Section',
}


# EXTENSION SETTINGS
myst_enable_checkboxes = True
myst_heading_anchors = 3
todo_include_todos = True

link_resolver_url = "https://github.com/thespacedoctor/soxspipe/blob/master"

remove_from_toctrees = ["utils/[!_]*"]


# BIBTEX STUFF
bibtex_bibfiles = ['dry-bookends-references.bib']
# bibtex_reference_style = 'author_year'
# bibtex_default_style = 'unsrtalpha'

# AUTODOC2 AND BIBTEX PLUGS CLASH DUE TO THIS ISSUE: https://github.com/pylint-dev/astroid/issues/2191
# UNTIL THIS IS FIXED SWITCH BETWEEN THE TWO
runAutodoc2 = os.getenv("AUTODOC2")


if True or (runAutodoc2 and runAutodoc2 != "None"):
    print("RUNNING AUTODOC2")
    autodoc2_packages = [
        {
            "path": "../../soxspipe",
            "exclude_files": ["*test_*.py"],

        }
    ]
    autodoc2_render_plugin = "myst"
    autodoc2_skip_module_regexes = [r".*test.*"]
    autodoc2_hidden_objects = ["private"]
    autodoc2_sort_names = True
# else:

if False:
    # THIS BUG NEEDS FIXED TO BE ABLE TO USE AUTODOC2 AND ROUND BRACKETS
    # https://github.com/pylint-dev/astroid/issues/2191
    print("RUNNING A SPHINX BUILD")
    import dataclasses
    from sphinxcontrib.bibtex.style.referencing.author_year import AuthorYearReferenceStyle
    from sphinxcontrib.bibtex.style.referencing import BracketStyle
    import sphinxcontrib.bibtex.plugin

    def bracket_style() -> BracketStyle:
        return BracketStyle(
            left='(',
            right=')',
        )

    @dataclasses.dataclass
    class MyReferenceStyle(AuthorYearReferenceStyle):
        bracket_parenthetical: BracketStyle = dataclasses.field(default_factory=bracket_style)
        bracket_textual: BracketStyle = dataclasses.field(default_factory=bracket_style)
        bracket_author: BracketStyle = dataclasses.field(default_factory=bracket_style)
        bracket_label: BracketStyle = dataclasses.field(default_factory=bracket_style)
        bracket_year: BracketStyle = dataclasses.field(default_factory=bracket_style)

    sphinxcontrib.bibtex.plugin.register_plugin(
        'sphinxcontrib.bibtex.style.referencing',
        'author_year_round', MyReferenceStyle)
    bibtex_reference_style = 'author_year_round'


# use the tab-trigger below for new function
# xt-def-function


# OpenGraph metadata
ogp_site_url = "https://soxspipe.readthedocs.io/en"
# This is the image that GitHub stores for our social media previews
ogp_image = "https://live.staticflickr.com/65535/51602359158_6105a9d0c7_b.png"
ogp_custom_meta_tags = [
    '<meta name="twitter:card" content="summary_large_image">',
]

# -- THEME SETTINGS -------------------------------------------------
html_theme = 'furo'
html_static_path = ['_static', '_images']
html_theme_options = {
    "light_logo": "thespacedoctor_icon_dark_circle.png",
    "dark_logo": "thespacedoctor_icon_white_circle.png",
    "source_repository": "https://github.com/thespacedoctor/soxspipe/",
    "source_branch": "master",
    "source_directory": "docs/source/",
}
html_favicon = "_images/favicon.ico"
html_title = f"soxspipe <small>v{__version__}</small>"
html_show_sourcelink = True
html_add_permalinks = u"  ∞"
# OTHER USEFUL SETTINGS
# html_theme_options = {
#     "announcement": "<em>Important</em> announcement!",
# }

# -- LaTeX output -------------------------------------------------

latex_engine = 'xelatex'

if not on_rtd:
    latex_documents = [
        ('overleaf/introduction', 'introduction.tex', u'introduction',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/recipes', 'recipes.tex', u'recipes',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/utils', 'utils.tex', u'utils',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/instruments', 'instruments.tex', u'instruments',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/data_organiser_and_reducer', 'data_organiser_and_reducer.tex', u'data organiser and reducer',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/data_reduction_cascades', 'data_reduction_cascades.tex', u'data reduction cascades',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/observing_modes', 'observing_modes.tex', u'observing modes',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/files', 'files.tex', u'files',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/quickstart_guide', 'quickstart_guide.tex', u'quickstart guide',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/installation', 'installation.tex', u'installation',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/preparing_a_workspace', 'preparing_a_workspace.tex', u'preparing a workspace',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/reductions', 'reductions.tex', u'reductions',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/pipeline_settings', 'pipeline_settings.tex', u'pipeline settings',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/logging', 'logging.tex', u'logging',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/sessions', 'sessions.tex', u'sessions',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/support', 'support.tex', u'sessions',
         u'David R. Young & Marco Landoni', 'howto', False),
        ('overleaf/references', 'references.tex', u'references',
         u'David R. Young & Marco Landoni', 'howto', False),
    ]
    latex_documents = [
        ('overleaf/recipes', 'recipes.tex', u'recipes',
         u'David R. Young & Marco Landoni', 'howto', False),
    ]
else:
    latex_documents = [
        ("index", 'soxspipe.tex', 'soxspipe Documentation', u'David R. Young & Marco Landoni', 'manual'),
    ]


latex_toplevel_sectioning = "section"
latex_show_urls = 'footnote'


def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None

    filename = info['module'].replace('.', '/')
    if info['fullname'] and "." not in info['fullname']:
        filename += "/" + info['fullname'] + ".py"
    else:
        if "/" in filename:
            filename = ("/").join(filename.split("/")[0:-1]) + "/"
        else:
            filename = ""
        filename += ("/").join(info['fullname'].split(
            ".")[0:-1]) + ".py" + "#" + info['fullname'].split(".")[-1]
    return link_resolver_url + "/" + filename


def updateUsageMd():
    """
    *Grab the usage from cl_utils.py to display in README.md*
    """
    from soxspipe import cl_utils
    import codecs
    usage = cl_utils.__doc__

    if not "Usage:" in usage or "todo:" in usage:
        return None
    usageString = ""
    for l in usage.split("\n"):
        usageString += "    " + l + "\n"

    usage = """

```bash
%(usageString)s
```
""" % locals()

    moduleDirectory = os.path.dirname(__file__)
    uFile = moduleDirectory + "/usage.md"
    writeFile = codecs.open(uFile, encoding='utf-8', mode='w')
    writeFile.write(usage)
    writeFile.close()

    return None


updateUsageMd()
