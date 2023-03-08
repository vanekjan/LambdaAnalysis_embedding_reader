#!/bin/csh

rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`

mkdir ./production/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}
mkdir ./jobs/report/${prodId}
mkdir ./jobs/list/${prodId}
mkdir ./jobs/csh/${prodId}
mkdir ./jobs/submit/${prodId}

 
star-submit-template -template submit.xml -entities productionId=${prodId},mFileList=${1},localPath=$PWD

mv *.session.xml ./jobs/submit/${prodId}
mv *.dataset ./jobs/submit/${prodId}
