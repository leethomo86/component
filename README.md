# component
Determine the coefficients that map two SCF determinants to each other.

Main/Component

NAME
      Determiant component calculator.
SYNOPSIS
      This program computes the components of a Slater determinant as an expansion (to first
      order) of another determinant.
USAGE
      berryCalc [-f1 <matrix_file_1>] [-f2 <matrix_file_2>] [--print-level <print_level>] [--method <method>]
        [--help]
OPTIONS
   1. Input/output

      -f1 matrix_file_1                First input file giving SCF wavefunction of determinant to
                                       be evaluated.

      -f2 matrix_file_2                Second input file giving SCF wavefunction of providing
                                       expansion for evaluation of matrix file 1.

      --print-level print_level        Verbosity of output. Default print level is 1. Options
                                       0-4.

      --sub-levels substitutions       Substitution levels permitted in truncated CI calculation.
                                       The default is all single substitutions.

                                       Example: [1,2] specifies single and double substitutions.

      --core-orbitals core-orbitals    Number of occupied orbitals to exclude from the truncated
                                       CI determinant expansion.

      --virt-orbitals virt-orbitals    Number of virtual orbitals to exclude from the truncated
                                       CI determinant expansion.

      --active-space active_space      Defines orbital space in which to construct initial basis
                                       determinant expansion through orbital swaps. For each input
                                       solution the number of electrons and orbitals in the active
                                       space should be specified.

                                       Example: [4,4:3,6] specifies an expansion of four electrons
                                       in four orbitals in the first input orbitals and an
                                       expansion of three electrons in six orbitals in the second
                                       input orbitals.

      --alter alter_list               Changes the order of orbitals in the input molecular
                                       orbitals. Alterations to each input orbitals should be
                                       colon separated and the two orbital numbers to be swapped
                                       should be separated by a if two alpha orbitals will be
                                       swapped, or b if two beta orbitals will be swapped.
                                       Different orbital swaps should be comma separated.

                                       Example: [3a4,2b4:5a6] swaps alpha orbitals 3 and 4 and
                                       beta orbitals 2 and 4 in input orbitals 1, and swaps alpha
                                       orbitals 5 and 6 in input orbitals 2.

      --ci-type type_string            Specifies the type of configuration interaction. Options
                                       are:
                                       1) oci (default)
                                          Perform orthogonal configuration interaction with
                                          determinant expansion specified by sub-levels option.
                                          Only one matrix file should be specified in the input
                                          file is expected and additional inputs will result in
                                          an error.
                                       2) ocas
                                          Perform orthogonal complete active space determinant
                                          expansion specified by active-space option. Only one
                                          matrix file should be specified in the input file is
                                          expected and additional inputs will result in an
                                          error.

      --help                           Output help documentation to terminal.

NOTES
      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
      f08 standard.

      Compilation tested using: gfortran 9.2.0

      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)

      Documentation generated with robodoc. To update documentation edit robodoc.rc to
      determine documentation output type and then run robodoc at the command line in the
      main directory.

    AUTHORS
      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
COPYRIGHT
      (c) 2022 by Lee M. Thompson distributed under terms of the MIT license.

---------------------------------------------------------------------------
