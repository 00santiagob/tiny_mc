# HP Pavillon

## Hardware

### CPU

```bash
❯ lscpu

Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         39 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  8
  On-line CPU(s) list:   0-7
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Core(TM) i5-9300H CPU @ 2.40GHz
    CPU family:          6
    Model:               158
    Thread(s) per core:  2
    Core(s) per socket:  4
    Socket(s):           1
    Stepping:            10
    CPU max MHz:         4100,0000
    CPU min MHz:         800,0000
    BogoMIPS:            4800.00
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mc
                         a cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss 
                         ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art
                          arch_perfmon pebs bts rep_good nopl xtopology nonstop_
                         tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cp
                         l vmx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1
                          sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsav
                         e avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_fault
                          epb invpcid_single pti ssbd ibrs ibpb stibp tpr_shadow
                          vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adjust 
                         bmi1 avx2 smep bmi2 erms invpcid mpx rdseed adx smap cl
                         flushopt intel_pt xsaveopt xsavec xgetbv1 xsaves dtherm
                          ida arat pln pts hwp hwp_notify hwp_act_window hwp_epp
                          md_clear flush_l1d arch_capabilities
Virtualisation features: 
  Virtualisation:        VT-x
Caches (sum of all):     
  L1d:                   128 KiB (4 instances)
  L1i:                   128 KiB (4 instances)
  L2:                    1 MiB (4 instances)
  L3:                    8 MiB (1 instance)
NUMA:                    
  NUMA node(s):          1
  NUMA node0 CPU(s):     0-7
Vulnerabilities:         
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushe
                         s, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Mitigation; IBRS
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer
                          sanitization
  Spectre v2:            Mitigation; IBRS, IBPB conditional, RSB filling, PBRSB-
                         eIBRS Not affected
  Srbds:                 Mitigation; Microcode
  Tsx async abort:       Not affected
```

### Periferics: Graphics

```
❯ lspci | grep -i vga

```

## Software

### Compilers

```bash
❯ dpkg -l | grep -w "compiler" | grep -w "C" | awk '/^ii/ {print $2}'

```