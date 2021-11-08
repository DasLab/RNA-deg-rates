function c = get_color_for_master_plot(RNA_type);
c = 'k';
if ~isempty(strfind(RNA_type, 'Roll your own structure 2' )); c = [1 0.5 0.5]; return; end;
if ~isempty(strfind(RNA_type, 'm1PSU' )); c = [0.1 0.1 0.7]; return; end;
if ~isempty(strfind(RNA_type, 'full mRNA PSU' )); c = [0.5 0.5 1]; return; end;
if ~isempty(strfind(RNA_type, 'combo18-PSU' )); c = [0.7 0.7 1]; return; end;
if ~isempty(strfind(RNA_type, 'combo18' )); c = [0.7 0.7 0.7]; return; end;
if ~isempty(strfind(RNA_type, 'PSU' )); c = [0.2 0.2 1]; return; end;
if ~isempty(strfind(RNA_type, 'm5C' )); c = [0.2 0.2 0.5]; return; end;
if ~isempty(strfind(RNA_type, '5mC' )); c = [0.2 0.2 0.5]; return; end;

if ~isempty(strfind(RNA_type, 'NLuc Barna' )); c = [0.8 0.1 0.1]; return; end;
if ~isempty(strfind(RNA_type, 'NLuc' )); c = [0.2 0.7 0.2]; end;
if ~isempty(strfind(RNA_type, 'Roll your own structure Followup' )); c = [0.7 0.5 0.1]; return; end;
if ~isempty(strfind(RNA_type, 'MEV' )); c = [1 0.3 0.8]; return; end;
if ~isempty(strfind(RNA_type, 'GFP' )); c = [1 0.5 0.3]; return; end;
if ~isempty(strfind(RNA_type, 'PSU' )); c = [1 0.5 0.3]; return; end;

if ~isempty(strfind(RNA_type, 'gel degradation' )); c = [0.3 0.3 0.3]; return; end;


