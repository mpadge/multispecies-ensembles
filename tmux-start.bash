#!/bin/sh
SESSION="four-node-nets"

tmux -2 new-session -d -s $SESSION

cd ./src/
tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim pop-fns.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe pop-fns.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model.c++' C-m

tmux split-window -h
tmux send-keys -t $SESSION:1 'vim model050-theoretical.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model050-theoretical.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model050.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model050.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model-envcor.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe model-envcor.c++' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -n utils
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim utils.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe utils.c++' C-m
tmux split-window -h
tmux send-keys -t $SESSION:2 'vim trophic-levels.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe trophic-levels.c++' C-m
tmux select-pane -t 0

cd ../
tmux new-window -t $SESSION:3 -n makefile
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim makefile' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe tmux-start.bash' C-m

tmux split-window -h
tmux select-pane -t 0

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
