"""Read-only localhost viewer for exactly one immutable delivery and measurement."""
import argparse,hashlib,json,sys,webbrowser
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))


def pick_entity_path():
    if sys.platform=='win32':
        import ctypes
        from ctypes import wintypes
        class OpenFile(ctypes.Structure):
            _fields_=[('lStructSize',wintypes.DWORD),('hwndOwner',wintypes.HWND),('hInstance',wintypes.HINSTANCE),
                ('lpstrFilter',wintypes.LPCWSTR),('lpstrCustomFilter',wintypes.LPWSTR),('nMaxCustFilter',wintypes.DWORD),
                ('nFilterIndex',wintypes.DWORD),('lpstrFile',wintypes.LPWSTR),('nMaxFile',wintypes.DWORD),
                ('lpstrFileTitle',wintypes.LPWSTR),('nMaxFileTitle',wintypes.DWORD),('lpstrInitialDir',wintypes.LPCWSTR),
                ('lpstrTitle',wintypes.LPCWSTR),('Flags',wintypes.DWORD),('nFileOffset',wintypes.WORD),
                ('nFileExtension',wintypes.WORD),('lpstrDefExt',wintypes.LPCWSTR),('lCustData',wintypes.LPARAM),
                ('lpfnHook',ctypes.c_void_p),('lpTemplateName',wintypes.LPCWSTR),('pvReserved',ctypes.c_void_p),
                ('dwReserved',wintypes.DWORD),('FlagsEx',wintypes.DWORD)]
        buffer=ctypes.create_unicode_buffer(32768);dialog=OpenFile();dialog.lStructSize=ctypes.sizeof(dialog)
        dialog.lpstrFile=ctypes.cast(buffer,wintypes.LPWSTR);dialog.nMaxFile=len(buffer)
        dialog.lpstrFilter='实体或三角网格 (*.json;*.ply)\0*.json;*.ply\0所有文件 (*.*)\0*.*\0\0';dialog.nFilterIndex=1
        dialog.lpstrInitialDir=str(Path.cwd())
        dialog.lpstrTitle='选择实体 JSON 或独立三角网格 PLY';dialog.Flags=0x00001000|0x00000800|0x00080000|0x00000008
        if ctypes.windll.comdlg32.GetOpenFileNameW(ctypes.byref(dialog)):return buffer.value
        error=ctypes.windll.comdlg32.CommDlgExtendedError()
        if error:raise RuntimeError('NATIVE_FILE_DIALOG_ERROR_'+str(error))
        return None
    import tkinter as tk
    from tkinter import filedialog
    root=tk.Tk();root.withdraw()
    try:return filedialog.askopenfilename(title='选择实体 JSON 或三角网格 PLY',filetypes=[('实体及网格','*.json *.ply')]) or None
    finally:root.destroy()

def serve_native(prefix, options, data):
    from root_model_delivery_exports import file_reference
    from root_model_delivery_validation import validation_for
    from root_model_delivery_downloads import display_validation,export_download
    from root_model_delivery_reconstruction import reconstruction_comparison,reconstruction_csv
    root=next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    def bound(ref):
        path=(prefix.parent/ref['path']).resolve()
        if not path.is_relative_to(root):raise ValueError('NATIVE_DELIVERY_FILE_OUTSIDE_EXPERIMENT')
        current=file_reference(path)
        if current['sha256']!=ref['sha256'] or current['bytes']!=ref['bytes']:
            raise ValueError('NATIVE_DELIVERY_FILE_CHANGED')
        return path
    routes={
        '/':(prefix.with_name(prefix.name+'_viewer.html'),'text/html; charset=utf-8'),
        '/delivery.json':(prefix.with_name(prefix.name+'_delivery.json'),'application/json'),
        '/measurement.json':(bound(data['measurement']),'application/json'),
        '/model_delivery_viewer.js':(P/'src/model_delivery_viewer.js','text/javascript; charset=utf-8'),
        '/graph.json':(bound(data['exports']['graph']),'application/json'),
        '/surface_positions.bin':(bound(data['surface']['buffers']['positions']),'application/octet-stream'),
        '/surface_indices.bin':(bound(data['surface']['buffers']['indices']),'application/octet-stream')}
    for name,ref in data['exports'].items():
        routes['/export/'+name]=(bound(ref),'application/json' if name=='graph' else 'text/csv; charset=utf-8')
    class NativeReadOnly(BaseHTTPRequestHandler):
        def do_GET(self):
            key=self.path.split('?')[0]
            if key in ['/reconstruction_comparison.json','/export/reconstruction']:
                try:
                    if key=='/reconstruction_comparison.json':
                        raw=json.dumps(reconstruction_comparison(prefix),ensure_ascii=False).encode('utf-8');mime='application/json'
                    else:raw=reconstruction_csv(prefix);mime='text/csv; charset=utf-8'
                except (ValueError,OSError,KeyError):self.send_error(409);return
                self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(raw)))
                if key.startswith('/export/'):self.send_header('Content-Disposition','attachment; filename="reconstruction_comparison.csv"')
                self.end_headers();self.wfile.write(raw);return
            if key in ['/parameters.json','/export/parameters','/export/parameters_json']:
                from root_model_delivery_downloads import parameter_download
                try:raw,mime,filename=parameter_download(prefix,'csv' if key=='/export/parameters' else 'json')
                except (ValueError,OSError,KeyError):self.send_error(409);return
                self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(raw)))
                if key.startswith('/export/'):self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
                self.end_headers();self.wfile.write(raw);return
            if key=='/validation.json':
                try:raw=json.dumps(display_validation(prefix),ensure_ascii=False).encode('utf-8')
                except (ValueError,OSError,KeyError):self.send_error(409);return
                self.send_response(200);self.send_header('Content-Type','application/json')
                self.send_header('Content-Length',str(len(raw)));self.end_headers();self.wfile.write(raw);return
            entry=routes.get(key)
            if entry is None:self.send_error(404);return
            if key in ['/export/paths','/export/samples','/export/graph','/graph.json']:
                kind=key.split('/')[-1] if key!='/graph.json' else 'graph'
                try:raw,mime,filename=export_download(prefix,kind)
                except (ValueError,OSError,KeyError):self.send_error(409);return
                self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(raw)))
                if key.startswith('/export/'):self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
                self.end_headers();self.wfile.write(raw);return
            path,kind=entry
            self.send_response(200);self.send_header('Content-Type',kind)
            self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data:; object-src 'none'; base-uri 'self'")
            self.send_header('Content-Length',str(path.stat().st_size))
            self.send_header('Cache-Control','no-store')
            if key.startswith('/export/'):
                self.send_header('Content-Disposition','attachment; filename="'+path.name+'"')
            self.end_headers()
            try:
                with path.open('rb') as stream:
                    for block in iter(lambda:stream.read(1048576),b''):self.wfile.write(block)
            except (BrokenPipeError,ConnectionResetError,ConnectionAbortedError):pass
        def do_POST(self):self.send_error(405)
        def log_message(self,*_):pass
    server=ThreadingHTTPServer(('127.0.0.1',options.port),NativeReadOnly)
    url='http://127.0.0.1:%d/'%server.server_address[1]
    print(url,flush=True)
    if options.open_browser:webbrowser.open(url)
    server.serve_forever()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prefix',type=Path);mode.add_argument('--workspace',type=Path)
    mode.add_argument('--pick-entity',action='store_true',help=argparse.SUPPRESS)
    ap.add_argument('--port',type=int,default=0)
    ap.add_argument('--open-browser',action='store_true')
    a=ap.parse_args()
    if a.pick_entity:
        print(json.dumps(dict(path=pick_entity_path()),ensure_ascii=False));return
    if a.workspace is not None:
        from root_model_delivery_service import serve_workspace
        return serve_workspace(a.workspace,a.port,a.open_browser)
    prefix=a.prefix.resolve()
    delivery=prefix.with_name(prefix.name+'_delivery.json');viewer=prefix.with_name(prefix.name+'_viewer.html');d=json.loads(delivery.read_text(encoding='utf-8-sig'))
    if d.get('schema_version')=='root-native-delivery-v001':
        return serve_native(prefix,a,d)
    root=next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file());measurement=(root/d['measurement']['experiment_relative_path']).resolve()
    assert measurement.is_relative_to(root) and hashlib.sha256(measurement.read_bytes()).hexdigest()==d['measurement']['sha256'] and measurement.stat().st_size==d['measurement']['bytes']
    routes={'/':(viewer,'text/html; charset=utf-8'),'/viewer.html':(viewer,'text/html; charset=utf-8'),'/delivery.json':(delivery,'application/json'),'/measurement.json':(measurement,'application/json')}
    model_preview=dict(available=False,status='ORIGINAL_MODEL_PREVIEW_NOT_BOUND')
    if d.get('original_model'):
        ref=d['original_model'];model=(delivery.parent/ref['relative_to_delivery']).resolve()
        assert model.stat().st_size==ref['bytes'] and hashlib.sha256(model.read_bytes()).hexdigest()==ref['sha256']
        if ref['bytes']<=ref['preview_byte_guard']:
            from root_triangle_measurement_input import read_triangle_ply_bytes
            try:
                mesh=read_triangle_ply_bytes(model.read_bytes(),ref['sha256'],max_vertices=3000,max_faces=1000)
                model_preview=dict(available=True,source_sha256=ref['sha256'],vertices=mesh.vertices.tolist(),faces=mesh.faces.tolist(),original_coordinates_unchanged=True,original_arrays_copied_to_file=False)
            except ValueError as exc:model_preview=dict(available=False,status='ORIGINAL_MODEL_PREVIEW_GUARDED',reason=str(exc))
        else:model_preview=dict(available=False,status='ORIGINAL_MODEL_PREVIEW_BYTE_GUARD')
    class ReadOnly(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path.split('?')[0]=='/model_preview.json':
                raw=json.dumps(model_preview,allow_nan=False).encode();self.send_response(200);self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(raw)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(raw);return
            entry=routes.get(self.path.split('?')[0])
            if entry is None:self.send_error(404);return
            path,kind=entry;raw=path.read_bytes();self.send_response(200);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(raw)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(raw)
        def do_POST(self):self.send_error(405)
        def log_message(self,*_):pass
    server=ThreadingHTTPServer(('127.0.0.1',a.port),ReadOnly);print('http://127.0.0.1:%d/'%server.server_address[1],flush=True);server.serve_forever()
if __name__=='__main__':main()
