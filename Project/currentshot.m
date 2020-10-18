function shotnum=currentshot()
mdsconnect('192.168.20.11');
mdsvalue('tcl($)','set tree EXL50');
shotnum=mdsvalue('current_shot(''EXL50'')');
mdsclose();
end