"""Drive the Streamlit script through all five steps with a stub `streamlit` module."""
import os, sys, types, runpy, shutil
from pathlib import Path
from unittest.mock import MagicMock
sys.path[:0] = [str(Path(__file__).resolve().parents[1]), str(Path(__file__).resolve().parent)]
from seqguide import demo as simulate
from seqguide.seqio import revcomp

data = Path('/tmp/sg_data'); shutil.rmtree(data, ignore_errors=True); (data/'input').mkdir(parents=True)
seq, feats = simulate.build_plasmid()
simulate.write_genbank(data/'input'/'pTest.gb', 'pTest', seq, feats)
simulate.sim_illumina_pairs(seq, 2000, data/'input'/'S1')
mut, _ = simulate.mutate(seq, feats); Path('/tmp/asm.txt').write_text(revcomp(mut))
os.environ.update(SEQGUIDE_DATA=str(data), FAKE_ASSEMBLY='/tmp/asm.txt',
                  PATH=str(Path(__file__).resolve().parent / 'fake_tools') + ':' + os.environ['PATH'])

class Rerun(Exception): pass
class Stop(Exception): pass
class SS(dict):
    __getattr__ = lambda s, k: s[k] if k in s else (_ for _ in ()).throw(AttributeError(k))
    __setattr__ = dict.__setitem__

st = types.ModuleType('streamlit'); st.session_state = SS()
def anyfn(*a, **k): return MagicMock()
for n in 'set_page_config title caption write divider header subheader markdown info success warning error metric image code progress empty spinner expander download_button html slider'.split():
    setattr(st, n, anyfn)
st.expander = lambda *a, **k: MagicMock()
st.spinner = lambda *a, **k: MagicMock()
st.slider = lambda label, lo, hi, val, **k: val
st.sidebar = MagicMock()
st.button = lambda label, **k: label.startswith(('Inspect', 'Continue', 'Start analysis')) and not k.get('disabled')
st.checkbox = lambda label, **k: label.endswith('MB)')
st.radio = lambda label, opts, **k: opts[0]
st.selectbox = lambda label, opts, **k: list(opts)[0]
st.number_input = lambda label, **k: k.get('value', 0)
st.file_uploader = lambda *a, **k: None
def _col():
    m = MagicMock(); m.button = st.button; return m
st.columns = lambda n: [_col() for _ in range(n if isinstance(n, int) else len(n))]
st.rerun = lambda: (_ for _ in ()).throw(Rerun())
st.stop = lambda: (_ for _ in ()).throw(Stop())
comp = types.ModuleType('streamlit.components'); v1 = types.ModuleType('streamlit.components.v1')
v1.html = lambda *a, **k: print('   embedded report:', len(a[0]), 'chars')
sys.modules.update({'streamlit': st, 'streamlit.components': comp, 'streamlit.components.v1': v1})

APP = str(Path(__file__).resolve().parents[1] / 'app' / 'streamlit_app.py')
for _ in range(80):
    step = st.session_state.get('step', 0)
    print('step', step)
    try:
        runpy.run_path(APP, run_name='__main__')
        if step == 4: break
    except Rerun:
        pass
    except Stop:
        print('   stopped'); break
print('final step', st.session_state['step'], 'error:', st.session_state['results']['error'])
print('variants:', [v['description'] for v in st.session_state['results']['comparison']['variants']])


# ---- second pass: browser mode (stlite) with no external tools -> built-in engine
import importlib
st.session_state.clear()
os.environ['PATH'] = '/usr/bin:/bin'
real_platform = sys.platform
for _ in range(80):
    step = st.session_state.get('step', 0)
    sys.platform = 'emscripten'
    try:
        runpy.run_path(APP, run_name='__main__')
        if step == 4: break
    except Rerun:
        pass
    except Stop:
        print('   stopped'); break
    finally:
        sys.platform = real_platform
r = st.session_state['results']
print('BROWSER MODE final step', st.session_state['step'], 'engine', r['params'].get('engine'), 'error:', r['error'])
print('variants:', [v['description'] for v in r['comparison']['variants']])
