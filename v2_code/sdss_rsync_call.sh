#rsync -avz --dry-run --include="*/" --include='frame-*.fits.bz2' --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ frames/
rsync -avz --include="*/" --include='frame-*.fits.bz2' --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ frames/
