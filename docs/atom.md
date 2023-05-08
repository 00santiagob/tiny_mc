# ASUS ZenBook Ux4

## Hardware

### CPU

```bash
❯ lscpu

Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         48 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  96
  On-line CPU(s) list:   0-47
  Off-line CPU(s) list:  48-95
Vendor ID:               AuthenticAMD
  Model name:            AMD EPYC 7643 48-Core Processor
    CPU family:          25
    Model:               1
    Thread(s) per core:  1
    Core(s) per socket:  48
    Socket(s):           1
    Stepping:            1
    Frequency boost:     enabled
    CPU(s) scaling MHz:  46%
    CPU max MHz:         3640.9170
    CPU min MHz:         0.0000
    BogoMIPS:            4600.42
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx mmxext fxsr_opt pdpe1gb rdtscp lm constant_tsc rep_good nopl nonstop_tsc cpuid extd_apicid aperfmperf rapl pni pclmulqdq monitor ssse3 fma cx16 pcid sse4_1 sse4_2 movbe popcnt aes xsave avx f16c rdrand lahf_lm cmp_legacy svm extapic cr8_legacy abm sse4a misalignsse 3dnowprefetch osvw ibs skinit wdt tce topoext perfctr_core perfctr_nb bpext perfctr_llc mwaitx cpb cat_l3 cdp_l3 invpcid_single hw_pstate ssbd mba ibrs ibpb stibp vmmcall fsgsbase bmi1 avx2 smep bmi2 invpcid cqm rdt_a rdseed adx smap clflushopt clwb sha_ni xsaveopt xsavec xgetbv1 xsaves cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local clzero irperf xsaveerptr rdpru wbnoinvd amd_ppin brs arat npt lbrv svm_lock nrip_save tsc_scale vmcb_clean flushbyasid decodeassists pausefilter pfthreshold v_vmsave_vmload vgif v_spec_ctrl umip pku ospke vaes vpclmulqdq rdpid overflow_recov succor smca
Virtualization features: 
  Virtualization:        AMD-V
Caches (sum of all):     
  L1d:                   1.5 MiB (48 instances)
  L1i:                   1.5 MiB (48 instances)
  L2:                    24 MiB (48 instances)
  L3:                    256 MiB (8 instances)
NUMA:                    
  NUMA node(s):          1
  NUMA node0 CPU(s):     0-47
Vulnerabilities:         
  Itlb multihit:         Not affected
  L1tf:                  Not affected
  Mds:                   Not affected
  Meltdown:              Not affected
  Mmio stale data:       Not affected
  Retbleed:              Not affected
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; Retpolines, IBPB conditional, IBRS_FW, STIBP always-on, RSB filling, PBRSB-eI
                         BRS Not affected
  Srbds:                 Not affected
  Tsx async abort:       Not affected
```

### Periferics: Graphics

```
❯ lspci | grep -i vga

44:00.0 VGA compatible controller: ASPEED Technology, Inc. ASPEED Graphics Family (rev 41)
```

## Software

### Compilers

```bash
❯ dpkg -l | grep -w "compiler" | grep -w "C" | awk '/^ii/ {print $2}'

aocc-compiler-3.2.0
aocc-compiler-4.0.0
clang
clang-11
clang-13
clang-14
clang-15
g++
g++-10
g++-11
g++-12
gcc
gcc-10
gcc-11
gcc-12
gcc-12-aarch64-linux-gnu
gcc-12-multilib
intel-oneapi-compiler-cpp-eclipse-cfg
intel-oneapi-compiler-dpcpp-cpp
intel-oneapi-compiler-dpcpp-cpp-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-common-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-classic-fortran-shared-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-common-2023.0.0
intel-oneapi-compiler-dpcpp-cpp-runtime-2023.0.0
intel-oneapi-compiler-dpcpp-eclipse-cfg
intel-oneapi-icc-eclipse-plugin-cpp-2023.0.0
```