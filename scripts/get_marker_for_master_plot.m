function s = get_marker_for_master_plot(RNA_type);
s = 'o';
if ~isempty(strfind(RNA_type, 'combo' )); s = '^'; return; end;
if ~isempty(strfind(RNA_type, 'final' )); s = 'd'; return; end;
if ~isempty(strfind(RNA_type, 'full' )); s = 's'; return; end;
