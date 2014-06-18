cd $PETSC_TMP
args *.log*

let g:logfilterpath = expand('<sfile>:p')

nnoremap <leader>r :execute ':so ' . g:logfilterpath<cr>

function! ErrFilter()
  g/PETSC ERROR: ---------/d
  g/PETSC ERROR: See docs/d
  g/PETSC ERROR: .* is given/d
  g/PETSC ERROR: Try option/d
  g/PETSC ERROR: or see/d
  g/PETSC ERROR: Note: The EXACT/d
  g/PETSC ERROR: Petsc Release Version/d
  g/PETSC ERROR: Configure run at/d
  g/PETSC ERROR: Configure options/d
  g/PETSC ERROR: Libraries linked from/d
  g/PETSC ERROR: likely location of problem/d
  g/PETSC ERROR: .* INSTEAD the line number of/d
  g/PETSC ERROR: .* on a .* named .* by .*/d
  g/PETSC ERROR: User provided function.* in unknown/d
endfunction

function! InfoFilter()
  g/....PetscCommDuplicate.*/d
  g/\v.{4}MatAssemblyEnd_SeqAIJ.*/d
endfunction


syn match logError /.*ERROR.*/
syn match DCell /.*FiberField.*/

hi def logError guifg=red
hi def DCell    guifg=green


