# DNA-Storage 
![logo2](https://github.com/Saarco99/DNA-Storage/assets/95081597/be750694-9ca7-4517-94ab-a87884028adf)
Primer Library Design for Random Access

Overview

This project focuses on designing a robust and efficient primer library for large-scale DNA storage systems. It builds upon the approach described in the article "Random Access in Large-Scale DNA Data Storage" published in Nature Biotechnology and aims to overcome the challenges related to biological constraints and scalability.
The project implements an optimized algorithm to generate and validate primers that meet strict criteria for DNA storage, including GC content, Hamming distance, absence of long homopolymers, sequence complementarity, and melting temperature. The solution integrates bioinformatics tools such as Primer3, NUPACK, and CD-HIT to enhance performance and ensure the robustness of the primer library.
Motivation

Data production is increasing exponentially, and traditional storage technologies face capacity, longevity, and environmental durability challenges. DNA storage offers a promising solution due to its high density, longevity, and energy efficiency. However, practical use is limited due to high costs and slow read/write times, and significant research is required to optimize DNA sequences used for storing digital information.
Goals

Primer Library Development: Develop an efficient algorithm to generate primers that meet biological and system constraints for random access DNA storage.
Efficiency: Improve the algorithm's performance to ensure it runs within a reasonable timeframe.
Maximization: Increase the number of usable primers stored within a DNA pool.
Algorithm Profiling: Analyze and profile the algorithm to understand its bottlenecks and scalability.
Approach

Initial Steps
Naive Algorithm: A basic Python implementation was tested to generate primers but proved too slow for practical use.
Bioinformatics Tools:
Primer3: Used for melting temperature calculation.
NUPACK: For secondary structure validation.
CD-HIT: For sequence similarity filtering.
Advanced Strategy
Primer Generation: N random primers are generated, followed by filtering against specific restrictions such as GC content, homopolymers, self-complementarity, and Hamming distance.
Intra-Primer and Inter-Primer Checks: The algorithm uses a trie-based system to efficiently filter primers based on Hamming distance, avoiding unnecessary comparisons.
Efficiency: By utilizing sets to store reverse complement patterns, the solution achieves faster inter-primer validation.
Key Biological Constraints
GC Content: Between 45% and 55%.
Long Homo-polymers: Disallowed to prevent reading errors.
Sequence Complementarity: Both intra-primer and inter-primer self-complementarity are minimized to avoid unintended interactions.
Hamming Distance: Ensures primers are sufficiently different to avoid misidentification.
Tools & Technologies

Python: Primary programming language for algorithm development.
Primer3: For melting temperature filtering.
NUPACK: For secondary structure filtering.
CD-HIT: For sequence similarity clustering.
Trie Data Structure: Used for efficient Hamming distance validation.
Results

This approach effectively scales for large datasets, providing a high number of validated primers suitable for random access DNA storage while maintaining biological integrity. The algorithm performs significantly faster compared to basic methods, especially in filtering primers for reverse complementarity and Hamming distance checks.
Future Work

Further work will focus on optimizing the runtime of the algorithm for larger datasets and exploring new methods for reducing primer errors during synthesis and sequencing.
Contributors

Saar Cohen
Shahar Trabelsi
Advisor: Dr. Sarel Cohen
Collaboration: Prof. Dalit Naor
References

Lee Organick, et al., "Random Access in Large-Scale DNA Data Storage", Nature Biotechnology, 2018.
Additional bioinformatics resources: DNA Storage Alliance


results:
![final_results_1and10mil.jpg](extra%2Ffinal_results_1and10mil.jpg)
![Filtration Pipeline.jpg](extra%2FFiltration%20Pipeline.jpg)
![sorted_unsorted.jpg](extra%2Fsorted_unsorted.jpg)
![timing.jpg](extra%2Ftiming.jpg)
![hamming_distance with trie-tree.jpg](extra%2Fhamming_distance%20with%20trie-tree.jpg)