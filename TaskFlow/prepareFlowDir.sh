if [ $1 = "msas" ]; then 
   PHOTONFLOWDIR=/home/mike/git_afterburner/AnalysisSoftware/TaskFlow
fi

echo $PHOTONFLOWDIR

mkdir Software
ln -sf $PHOTONFLOWDIR/*.txt .
ln -sf $PHOTONFLOWDIR/*.sh .
ln -sf $PHOTONFLOWDIR/Software/*.C Software/
ln -sf $PHOTONFLOWDIR/Software/*.h Software/
