" since env vars not sourced:
let g:petsc_dir="/home/abergman/apps/petsc"

set nocscopeverbose
cscope add cscope.out
cscope add ~/apps/openmpi-1.8.1/build/include/cscope.out ~/apps/openmpi-1.8.1/build/include/
cscope add ~/apps/petsc/cscope.out ~/apps/petsc
set cscopeverbose
set tags+=~/apps/openmpi-1.8.1/build/include/tags
set tags+=~/apps/petsc/CTAGS
set path+=~/apps/petsc
set path+=~/apps/openmpi-1.8.1

map <leader>0 :cd ~/projects/DCell<CR>
map <leader>1 :CtrlP ~/apps/petsc<CR>
map <leader>2 :find local.vimrc<CR>
map <leader>3 :find /home/abergman/tmp/info.log.0<CR>
map <leader>4 :find ./sims/Fibers/Fibers.c<CR>

augroup setDCellDir
  autocmd!
  autocmd InsertLeave * :cd ~/projects/DCell
augroup END

augroup setNomodPETSc
  autocmd!
  autocmd BufEnter ~/apps/petsc/* setlocal nomodifiable 
augroup END
"if exists("g:did_localvimrc")
"  finish
"endif
"let g:did_localvimrc = 1
" stuff to lead only once

"let &efm = "%m line %l in %f,%-G[0]PETSC ERROR:%m" . efm

