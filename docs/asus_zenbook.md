# ASUS ZenBook UX410U

## Hardware

### CPU

```bash
❯ lscpu

Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         39 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  4
  On-line CPU(s) list:   0-3
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Core(TM) i7-7500U CPU @ 2.70GHz
    CPU family:          6
    Model:               142
    Thread(s) per core:  2
    Core(s) per socket:  2
    Socket(s):           1
    Stepping:            9
    CPU max MHz:         3500,0000
    CPU min MHz:         400,0000
    BogoMIPS:            5799.77
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_fault epb invpcid_single pti ssbd ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adjust bmi1 avx2 smep bmi2 erms invpcid mpx rdseed adx smap clflushopt intel_pt xsaveopt xsavec xgetbv1 xsaves dtherm ida arat pln pts hwp hwp_notify hwp_act_window hwp_epp md_clear flush_l1d arch_capabilities
Virtualization features: 
  Virtualization:        VT-x
Caches (sum of all):     
  L1d:                   64 KiB (2 instances)
  L1i:                   64 KiB (2 instances)
  L2:                    512 KiB (2 instances)
  L3:                    4 MiB (1 instance)
NUMA:                    
  NUMA node(s):          1
  NUMA node0 CPU(s):     0-3
Vulnerabilities:         
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushes, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Mitigation; IBRS
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; IBRS, IBPB conditional, RSB filling, PBRSB-eIBRS Not affected
  Srbds:                 Mitigation; Microcode
  Tsx async abort:       Not affected
```

### Periferics: Graphics

```
❯ lspci | grep -i vga

00:02.0 VGA compatible controller: Intel Corporation HD Graphics 620 (rev 02)
```

## Software

### Compilers

```bash
❯ dpkg -l | grep -w "compiler" | grep -w "C" | awk '/^ii/ {print $2}'
clang
clang-14
g++
g++-11
gcc
gcc-11
intel-oneapi-compiler-cpp-eclipse-cfg
intel-oneapi-compiler-dpcpp-cpp
intel-oneapi-compiler-dpcpp-cpp-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-common-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-common-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-runtime-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-classic-fortran-shared-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-classic-fortran-shared-runtime-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-common-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-common-2023.1.0
intel-oneapi-compiler-dpcpp-cpp-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-runtime-2023.1.0
intel-oneapi-compiler-dpcpp-eclipse-cfg
intel-oneapi-icc-eclipse-plugin-cpp-2023.0.0
intel-oneapi-icc-eclipse-plugin-cpp-2023.1.0
tcc
```