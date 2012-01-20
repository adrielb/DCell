rsync -avh --progress --stats --delete --checksum --inplace \
--exclude 'levlabnas/' \
--exclude 'Music/' \
--exclude '.mozilla/firefox/gz686kye.default/Cache/' \
/home/abergman/ levlabnas:whisperquiet

