module test_module
    implicit none
contains
    subroutine test_mv()

        real(kind = 8) :: test_array(3,3)[*]

        test_array = 1.

