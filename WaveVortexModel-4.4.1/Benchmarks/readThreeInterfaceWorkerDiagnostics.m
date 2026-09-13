function diagnostics = readThreeInterfaceWorkerDiagnostics(stdoutPath,stderrPath)
% Retain bounded worker output in receipts; failed workloads keep full files.
diagnostics = struct("stdout",readExcerpt(stdoutPath),"stderr",readExcerpt(stderrPath));
end

function record = readExcerpt(pathname)
limitBytes = 64*1024;
record = struct("text","","originalBytes",0,"truncated",false);
if ~isfile(pathname), return, end
[fileId,message] = fopen(pathname,"r");
if fileId < 0
    error("WaveVortexBenchmark:WorkerLogReadFailed","Unable to read %s: %s",pathname,message);
end
cleanup = onCleanup(@()fclose(fileId));
fseek(fileId,0,"eof");
record.originalBytes = ftell(fileId);
fseek(fileId,0,"bof");
record.truncated = record.originalBytes > limitBytes;
if record.truncated
    head = fread(fileId,limitBytes/2,"*uint8").';
    fseek(fileId,-limitBytes/2,"eof");
    tail = fread(fileId,limitBytes/2,"*uint8").';
    record.text = string(native2unicode(head,"UTF-8"))+newline+"[worker log truncated; full file retained on failure]"+newline+string(native2unicode(tail,"UTF-8"));
else
    record.text = string(native2unicode(fread(fileId,Inf,"*uint8").',"UTF-8"));
end
end
