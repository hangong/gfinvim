try
    echo on
    mex dval.c
    mex fgt_model.c
    mex fgt_predict.c
    echo off
catch ME
    idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
    if(~isempty(idSegLast))
        disp('Mex compilation error, please be sure to setup your compiler by typing ''mex -setup'' ')
    end
end