      program component
!
!     This program computes the components of a Slater determinant as an expansion (to first 
!     order) of another determinant.
!
!     L. M. Thompson, 2022
!
      use MQC_Gaussian
      use omp_lib
      use iso_fortran_env, only: int32, int64, real64
!
!****x* Main/Component
!*    NAME
!*      Determiant component calculator. 
!*
!*    SYNOPSIS
!*      This program computes the components of a Slater determinant as an expansion (to first 
!*      order) of another determinant.
!
      implicit none
      character(len=:),allocatable::command,fileName1,fileName2,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo
      integer(kind=int64)::iOut=output_unit,iIn=input_unit,iPrint=1,iUnit,flag,stat_num,nullSize,&
        i,j,nCore=0,nVirt=0
      character(len=256)::sub_string='',activeSpace='',alterations='',ci_string='oci'
      type(mqc_scf_integral)::overlap,mo1
      type(mqc_scalar)::NIJ,pNIJ
      real(kind=real64),parameter::zero=0.0d0,half=0.5d0,one=1.0d0,zero_thresh=1.0d-8
      type(mqc_molecule_data)::molInfo
      type(mqc_scf_integral),dimension(:),allocatable::mo_list
      logical::UHF
      type(mqc_wavefunction)::wavefunction
      integer(kind=int64),dimension(:),allocatable::isubs,inactiveList,activeList,alphaList,betaList
      type(mqc_determinant)::determinants
      type(mqc_vector)::subs
!
!*    USAGE
!*      berryCalc [-f1 <matrix_file_1>] [-f2 <matrix_file_2>] [--print-level <print_level>] [--method <method>] 
!*        [--help] 
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  Determiant Component Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 22.06.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2022.',NEW_LINE('a')
!
!     Parse input options.
!
!*   1. Input/output
!*
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f1') then
!
!*      -f1 matrix_file_1                First input file giving SCF wavefunction of determinant to
!*                                       be evaluated.
!*
          call mqc_get_command_argument(i+1,fileName1)
          j = i+2
        elseif(command.eq.'-f2') then
!
!*      -f2 matrix_file_2                Second input file giving SCF wavefunction of providing 
!*                                       expansion for evaluation of matrix file 1.
!*
          call mqc_get_command_argument(i+1,fileName2)
          j = i+2
        elseif(command.eq.'--print-level') then
!
!*      --print-level print_level        Verbosity of output. Default print level is 1. Options
!*                                       0-4.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I1)') iPrint
          j = i + 2
        elseIf(command.eq.'--sub-levels') then
!
!*      --sub-levels substitutions       Substitution levels permitted in truncated CI calculation. 
!*                                       The default is all single substitutions.
!*
!*                                       Example: [1,2] specifies single and double substitutions.
!*
          call mqc_get_command_argument(i+1,command)
          sub_string = command
          j = i+2
        elseIf(command.eq.'--core-orbitals') then
!
!*      --core-orbitals core-orbitals    Number of occupied orbitals to exclude from the truncated 
!*                                       CI determinant expansion.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I3)') nCore
          j = i+2
        elseIf(command.eq.'--virt-orbitals') then
!
!*      --virt-orbitals virt-orbitals    Number of virtual orbitals to exclude from the truncated 
!*                                       CI determinant expansion.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I3)') nVirt
          j = i+2
        elseIf(command.eq.'--active-space') then
!
!*      --active-space active_space      Defines orbital space in which to construct initial basis
!*                                       determinant expansion through orbital swaps. For each input
!*                                       solution the number of electrons and orbitals in the active
!*                                       space should be specified.
!*
!*                                       Example: [4,4:3,6] specifies an expansion of four electrons
!*                                       in four orbitals in the first input orbitals and an
!*                                       expansion of three electrons in six orbitals in the second
!*                                       input orbitals.
!*
          call mqc_get_command_argument(i+1,command)
          activeSpace = command
          j = i+2
        elseIf(command.eq.'--alter') then
!
!*      --alter alter_list               Changes the order of orbitals in the input molecular
!*                                       orbitals. Alterations to each input orbitals should be
!*                                       colon separated and the two orbital numbers to be swapped
!*                                       should be separated by a if two alpha orbitals will be
!*                                       swapped, or b if two beta orbitals will be swapped.
!*                                       Different orbital swaps should be comma separated.
!*
!*                                       Example: [3a4,2b4:5a6] swaps alpha orbitals 3 and 4 and
!*                                       beta orbitals 2 and 4 in input orbitals 1, and swaps alpha
!*                                       orbitals 5 and 6 in input orbitals 2.
!*
          call mqc_get_command_argument(i+1,command)
          alterations = command
          j = i+2
        elseIf(command.eq.'--ci-type') then
!
!*      --ci-type type_string            Specifies the type of configuration interaction. Options
!*                                       are:
!*                                       1) oci (default)
!*                                          Perform orthogonal configuration interaction with 
!*                                          determinant expansion specified by sub-levels option.
!*                                          Only one matrix file should be specified in the input 
!*                                          file is expected and additional inputs will result in 
!*                                          an error.
!*                                       2) ocas
!*                                          Perform orthogonal complete active space determinant 
!*                                          expansion specified by active-space option. Only one 
!*                                          matrix file should be specified in the input file is 
!*                                          expected and additional inputs will result in an 
!*                                          error.
!*
          call mqc_get_command_argument(i+1,command)
          ci_string = command
          j = i+2
        elseIf(command.eq.'--help') then
!
!*      --help                           Output help documentation to terminal.
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count())
          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/component.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
      endDo
!
!     Parse input file and extract required data from matrix files.
!
      call fileInfo%getESTObj('mo coefficients',est_integral=mo1,filename=fileName1)
      call fileInfo%getESTObj('overlap',est_integral=overlap,filename=fileName1)
      allocate(mo_list(1))
      call fileInfo%getESTObj('mo coefficients',est_integral=mo_list(1),filename=fileName2)
      if(iPrint.ge.4) then
        call MO1%print(iOut,'MO coefficients from matrix file 1')
        call MO_list(1)%print(iOut,'MO coefficients from matrix file 2')
        call overlap%print(iOut,'Overlap from matrix file 1')
      endIf

      call fileInfo%getMolData(molInfo,filename=fileName1)
      call molInfo%print(iOut)
      call mqc_print(mqc_get_nuclear_repulsion(molInfo),iOut,'Vnn')

      call fileInfo%load(fileName2)
      call fileInfo%getESTObj('wavefunction',wavefunction)

      if(len_trim(alterations).ne.0) then
        call orbital_swapper(alterations,mo_list)
        if(iPrint.ge.4) call mo_list(1)%print(iOut,'Altered MO2 coefficients')
      endIf

      if(MQC_Gaussian_IsUnrestricted(fileInfo)) then
        UHF = .true.
        if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found UHF wavefunction'//NEW_LINE('A')
      elseIf(MQC_Gaussian_IsRestricted(fileInfo)) then
        UHF = .False.
        if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found RHF wavefunction'//NEW_LINE('A')
      else
        call mqc_error('General wavefunction type not supported for MO2')
      endIf 
      if (MQC_Gaussian_IsComplex(fileInfo)) call mqc_error('Complex wavefunctions unsupported')
!
!     Generate Slater determinants wavefunction expansion if orthogonal expansion requested. 
!
      if(ci_string.eq.'oci') then
        call substitution_builder(sub_string,subs)
        if(iPrint.ge.1) then
          write(iOut,'(1X,A)') 'Building Determinant Strings'
          call subs%print(6,'Permitted substitution levels',Blank_At_Bottom=.true.)
          write(iOut,'(1X,A,1X,I3,1X,A)') 'Excluding',nCore,'core orbitals'
          write(iOut,'(1X,A,1X,I3,1X,A)') 'Excluding',nVirt,'virtual orbitals'
        endIf
        if(nCore.gt.min(int(wavefunction%nAlpha),int(Wavefunction%nBeta))) &
          call mqc_error_i('Impossible number of core orbitals requested in truncated CI expansion',&
          6,'nCore',nCore,'wavefunction%nAlpha',int(wavefunction%nAlpha),'wavefunction%nBeta',&
          int(wavefunction%nBeta))
        if(nVirt.gt.min(int(wavefunction%nBasis-wavefunction%nAlpha),int(wavefunction%nBasis-Wavefunction%nBeta))) &
          call mqc_error_i('Impossible number of virtual orbitals excluded from truncated CI expansion',&
          6,'nVirt',nVirt,'wwavefunction%nBasis-avefunction%nAlpha',int(wavefunction%nBasis-wavefunction%nAlpha),&
          'wavefunction%nBasis-wavefunction%nBeta',int(wavefunction%nBasis-wavefunction%nBeta))
        isubs = [(i, i=1,int(maxval(subs)))]
        call trci_dets_string(iOut,iPrint,wavefunction%nBasis-nCore-nVirt,wavefunction%nAlpha-nCore, &
          Wavefunction%nBeta-nCore,isubs,determinants,nCore)
      elseIf(ci_string.eq.'ocas') then
        call parse_active_space(activeSpace,1,wavefunction%nBasis,wavefunction%nAlpha,&
          Wavefunction%nBeta,Wavefunction%nElectrons,activeList,inactiveList,alphaList,betaList)
        if(iPrint.ge.1) then
          write(iOut,'(1X,A)') 'Building Determinant Strings'
          call mqc_print(6,activeList,'Active orbitals')
          call mqc_print(6,inactiveList,'Core orbitals')
          call mqc_print(6,alphaList,'Active alpha electrons')
          call mqc_print(6,betaList,'Active beta electrons')
        endIf
        call gen_det_str(iOut,iPrint,activeList(1),alphaList(1),betaList(1),determinants,inactiveList(1))
        nCore = inactiveList(1)
        nVirt = wavefunction%nBasis-nCore-activeList(1)
      endIf

      do i = 1, determinants%NAlpStr
        do j = 1, determinants%NBetStr
          call get_nij(NIJ,pnIJ,nullSize,&
            mo_list(1)%orbitals('useStrings',MQC_Bit2Num_String(determinants%strings%alpha%vat([i],[0])),&
            MQC_Bit2Num_String(determinants%strings%beta%vat([j],[0]))),mo1,wavefunction%overlap_matrix,&
            wavefunction%nBasis,wavefunction%nAlpha,wavefunction%nBeta)
          call mqc_print(NIJ,iOut,&
            mqc_detString_print(determinants%strings%alpha%vat([i],[0]),determinants%strings%beta%vat([j],[0]),&
            int(wavefunction%nBasis),nCore,nVirt))
        endDo
      endDo
        
      
      contains

!
!     PROCEDURE get_nij
!     
!     get_nij is a subroutine that returns the overlap and psuedooverlap of two
!     (potentially nonorthogonal) Slater determinants, as well as the dimension
!     of the overlap null space, the overlap and psuedo-overlap matrix elements. The
!     sign of nullSize gives the multiple of matrix elements accounting for antisymmetry
!     due to permutation of orbitals in the SVD.
!
      subroutine get_nij(nIJ,pnIJ,nullSize,mo_I,mo_J,overlap,nBasis,nAlpha,nBeta)

      implicit none

!     input/output variables
      type(mqc_scalar),intent(inOut)::nIJ,pnIJ
      integer(kind=int64)::nullSize
      type(mqc_scf_integral),intent(in)::mo_I,mo_J,overlap
      type(mqc_scalar),intent(in)::nBasis,nAlpha,nBeta

!     subroutine variables
      type(mqc_matrix)::mo_I_occ,mo_J_occ,mIJ,uMat,vMat
      type(mqc_vector)::sigmaMat
      logical::orthflag
      integer(kind=int64)::i

!
      mo_I_occ = mqc_integral_output_block(mo_I%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      mo_J_occ = mqc_integral_output_block(mo_J%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      
      mIJ = matmul(matmul(dagger(mo_I_occ),overlap%getBlock('full')),mo_J_occ)
!      call mIJ%print(6,'LMTLMT mij')
      nIJ = mIJ%det()

      orthflag = .false.
      if((nIJ%abs()).lt.zero_thresh) then
        call mIJ%svd(EVals=sigmaMat,EUVecs=uMat,EVVecs=vMat)
        if(minval(abs(sigmaMat)).lt.zero_thresh) orthflag = .true.
      endIf

      if(orthflag) then
        pnIJ = 1.0
        nullSize = 0
        do i = 1,size(sigmaMat)
          if(sigmaMat%at(i).gt.zero_thresh) then
            pnIJ = pnIJ*sigmaMat%at(i)
          else
            nullSize = nullSize + 1
          endIf
        endDo
        nullSize = sign(1.0,real(uMat%det()))*sign(1.0,real(vMat%det()))*nullSize
      else
        pnIJ = nIJ
        nullSize = 0
      endIf
!
      end subroutine get_nij
!
!
      subroutine compute_NIJ(detLen,bra_weights,ket_weights,NIJ_matrix,NIJ)
!
      implicit none
      type(mqc_vector),intent(in)::bra_weights,ket_weights
      type(mqc_matrix),intent(in)::NIJ_matrix
      integer,intent(in)::detLen
      type(mqc_scalar),intent(out)::NIJ
      integer::i,j
      real(kind=real64),parameter::zero=0.0d0
!
      NIJ = zero
      do i = 1, max(detLen,1)
        do j = 1, max(detLen,1)
          NIJ = NIJ + conjg(bra_weights%at(i))*ket_weights%at(j)*NIJ_matrix%at(i,j)
        endDo
      endDo
!
      end subroutine compute_NIJ
!
!
!     PROCEDURE orbital_swapper
!
!     orbital_swapper is a subroutine that parses input and swaps molecular orbitals on the input matrix files.
!
      subroutine orbital_swapper(alterations,mo_list)

      implicit none

!     input/output variables
      type(mqc_scf_integral),dimension(:),allocatable,intent(inOut)::mo_list
      character(len=*),intent(in)::alterations

!     text parsing variables
      integer::i,leftNum,rightNum,solutionNum
      character(len=:),allocatable::alterString
      logical::lNum,rNum,lSet,rSet,betaPair,newSoln,hasNum
      character(len=80)::processString

      alterString = trim(alterations)
      if(len(alterString).ne.0) then
        !  Initialize orbital swap variables
        leftNum = 0
        rightNum = 0
        !  Flags for determining current state of the string processing.  lNum and rNum
        !  indicate whether the orbital pair is being assembled, lSet and rSet inidicate
        !  whether the left and right orbital numbers have been processed and stored.
        !  processString contains the characters that are being examined.
        lNum = .false.
        rNum = .false.
        lSet = .false.
        rSet = .false.
        betaPair = .false.
        newSoln = .false.
        hasNum = .false.
        processString = ''
        solutionNum = 1

        do i = 1, len(alterString)
          select case(alterString(i:i))
          case('[')
            if(i.eq.1) then
              lNum = .true.
              cycle
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            endif
          case("0":"9")
            if(i.lt.len(alterString)) then
              processString = trim(processString)//alterString(i:i)
              if(.not.hasNum) hasNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            endif
          case("a")
            if(lNum.and.hasNum) then
              read(processString,'(I3)') leftNum
              processString = ''
              betaPair = .false.
              hasNum = .false.
              lSet = .true.
              lNum = .false.
              rNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case("b")
            if(lNum.and.hasNum) then
              read(processString,'(I3)') leftNum
              processString = ''
              betaPair = .true.
              hasNum = .false.
              lSet = .true.
              lNum = .false.
              rNum=.true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case(",")
            if((lSet.and.rNum).and.hasNum) then
              read(processString,'(I3)') rightNum
              processString = ''
              hasNum = .false.
              rNum = .false.
              rSet = .true.
              lNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case(":")
            if((lSet.and.rNum).and.hasNum) then
              read(processString,'(I3)') rightNum
              processString = ''
              hasNum = .false.
              lNum = .true.
              rNum = .false.
              rSet = .true.
              newSoln = .true.
            else if((.not.lSet.and..not.rSet.and.lNum).and..not.(hasNum)) then
              newSoln = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case("]")
            if(rNum.and.hasNum) then
              read(processString,'(I3)') rightNum
              rNum = .false.
              rSet = .true.
            else if(.not.(alterString(i-1:i-1).eq.':')) then
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case default
            call mqc_error("Unrecognized character in alteration string input:  "//alterString(i:i))
          end select

          if(lSet.and.rSet) then
            if(.not.betaPair) then
              mo_list(solutionNum) = mo_list(solutionNum)%swap([leftNum,rightNum])
            else
              mo_list(solutionNum) = mo_list(solutionNum)%swap(betaOrbsIn=[leftNum,rightNum])
            endIf
            rSet = .false.
            lSet = .false.
          end if
          if(newSoln) then
            solutionNum = solutionNum + 1
            if(solutionNum.gt.size(mo_list)) call mqc_error_i('Orbital alterations requested on &
              &more determinants than available',6,'SolutionNum',SolutionNum,'size(mo_list)',&
              size(mo_list))
            newSoln = .false.
          end if
        end do
      end if

      end subroutine orbital_swapper
!
!
!     PROCEDURE substitution_builder
!
!     substitution_builder is a subroutine that builds the input wavefunction substituion levels.
!
      subroutine substitution_builder(subs_in,subs_out)

      implicit none

!     input/output variables
      character(len=*),intent(in)::subs_in
      type(mqc_vector),intent(inOut)::subs_out

!     text parsing variables
      character(len=:),allocatable::subsString
      character(len=10)::val
      integer::i,ival
      logical::newNum=.false.

      subsString = trim(subs_in)
      if(len(subsString).eq.0) then
        write(iOut,'(1X,A)') 'Defaulting to CIS.'//NEW_LINE('A')
        subs_out = [1]
      else
        do i = 1,len(subsString)
          select case (subsString(i:i))
          case('[','\(')
            if(i.eq.1) then
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            endIf
          case(']','\)')
            if(i.eq.len(subsString).and..not.newNum.and.i.ne.1) then
              read(val,'(I10)') ival
              call subs_out%push(ival)
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            endIf
          case(',',' ')
            if(i.eq.1.or.i.eq.len(subsString).or.i.eq.2.or.i.eq.len(subsString)-1.or.newNum) then
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            else
              read(val,'(I10)') ival
              call subs_out%push(ival)
              newNum = .true.
              cycle
            endIf
          case('0':'9')
            if(i.eq.1.or.i.eq.len(subsString)) then
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            else
              if(newNum) then
                val = subsString(i:i)
              else
                val = trim(val)//subsString(i:i)
              endIf
              newNum = .false.
              cycle
            endIf
          case default
            call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
          end select
        endDo
      endIf

      end subroutine substitution_builder
!    
!
!     PROCEDURE parse_active_space
!
!     parse_active_space is a subroutine that parses input for complete active space determinant 
!     expansions.
!
      subroutine parse_active_space(active_space,numFile,nBasis,nAlpha,nBeta,nElec,activeList,&
        inactiveList,alphaList,betaList)

      implicit none

!     input/output variables
      character(len=*),intent(in)::active_space
      integer,intent(in)::numFile
      type(mqc_scalar),intent(in)::nBasis,nElec,nAlpha,nBeta
      integer,dimension(:),allocatable,intent(out)::inactiveList,activeList,alphaList,betaList

!     text parsing variables
      integer::i,occNum,elecNum,solutionNum,maxActive
      character(len=:),allocatable::activeString
      logical::occSet,elecSet,oNum,eNum
      character(len=80)::processString
!
!     parse active space information
!
      allocate(inactiveList(numFile))
      allocate(activeList(numFile))
      allocate(alphaList(numFile))
      allocate(betaList(numFile))
      activeString = trim(active_space)
      if(len(activeString).eq.0) then
        do i = 1, numFile
          activeList(i) = nBasis
          inactiveList(i) = 0
          alphaList(i) = nAlpha
          betaList(i) = nBeta
          write(6,'(A)') ' Defaulting to substitutions over all orbitals in solution '//trim(num2char(i))
        endDo
      else
        !  Initialize occupied orbitals and electron count to zero
        occNum = 0
        elecNum = 0
        !  Flags for determining current state of the string processing.  oNum
        !  and eNum indicate whether the occupied orbital or active electron
        !  number for solution 'n' is currently being assembled.  occSet and
        !  elecSet indicate whether the occupied orbital number or active
        !  electron number have been fully assembled and stored.  processString
        !  contains the characters that correspond to occNum or elecNum.
        oNum = .false.
        eNum = .false.
        occSet = .false.
        elecSet = .false.
        processString = ''
        solutionNum = 1
        do i = 1, len(activeString)
          select case(activeString(i:i))
          case('[')
            if(i.eq.1) then
              oNum = .true.
              cycle
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            end if
          case("0":"9")
            if(i.lt.len(activeString)) then
              processString = trim(processString)//activeString(i:i)
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            endif
          case(",")
            if(.not.oNum) then
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(I3)') occNum
              processString = ''
              occSet = .true.
              oNum = .false.
              eNum = .true.
            end if
          case("]")
            if(.not.eNum) then
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(I3)') elecNum
              elecSet = .true.
              eNum = .false.
            end if
          case(":")
            if(occSet.and.eNum.and.(.not.elecSet)) then
              read(processString,'(I3)') elecNum
              elecSet = .true.
              oNum = .true.
              eNum = .false.
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            end if
          case default
            call mqc_error("Unrecognized character in active space input:  "//activeString(i:i))
          end select

          !  Algorithm has identified a valid pair of numbers for processing.
          if(occSet.and.elecSet) then

            maxActive = nBasis - (nElec/2) - mod(int(nElec),2) + (elecNum / 2) + mod(elecNum,2)

            !  Check users requested orbitals or electrons are sensible
            !  actually exist in the SCF results. I would hope that's not something
            !  the user would forget to check, but I'll err on the side of caution
            if(elecNum.gt.nElec) then
              call mqc_error("User has requested more active electrons than exist for solution "//&
                num2char(solutionNum))
            else if(occNum.gt.nBasis) then
              call mqc_error("User has requested more active orbitals than exist for solution "//&
                num2char(solutionNum))
            else if(occNum.gt.maxActive) then
              call mqc_error("User has requested too many active orbitals for solution "//&
                num2char(solutionNum))
            end if
   
            activeList(solutionNum) = occNum
            inactiveList(solutionNum) = (nelec-elecNum)/2 + mod(int(nelec-elecNum),2)
            alphaList(solutionNum) = nAlpha - inactiveList(solutionNum)
            betaList(solutionNum) = nBeta - inactiveList(solutionNum)
   
            occSet = .false.
            elecSet = .false.
            oNum = .true.
            eNum = .false.
            processString = ''
            solutionNum = solutionNum + 1
   
          end if
          if((solutionNum-1).gt.numFile) call mqc_error_i('Active space specifies more input solutions than&
            & present',6,'solutionNum',solutionNum-1,'numFile',numFile)
        end do
        if((solutionNum-1).lt.numFile) call mqc_error_i('Active space specifies fewer input solutions than&
          & present',6,'solutionNum',solutionNum-1,'numFile',numFile)
      end if

      end subroutine parse_active_space
!    
!    
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
!*      f08 standard.
!*
!*      Compilation tested using: gfortran 9.2.0
!*
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to
!*      determine documentation output type and then run robodoc at the command line in the
!*      main directory.
!*
!*    AUTHORS
!*      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2022 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
      end program component
