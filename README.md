[Primer Library Design for Random Access DNA Storage
](src/benchmarks/final_benchmark_full_pipeline.py )
![logo2](https://github.com/Saarco99/DNA-Storage/assets/95081597/be750694-9ca7-4517-94ab-a87884028adf)

### **Overview**

This project focuses on designing a robust and efficient primer library for large-scale DNA storage systems. It builds upon the approach described in the article "Random Access in Large-Scale DNA Data Storage" published in Nature Biotechnology and aims to overcome the challenges related to biological constraints and scalability.

#### **Goals**


Primer Library Development: Develop an efficient algorithm to generate primers that meet biological and system constraints for random access DNA storage.
Efficiency: Improve the algorithm's performance to ensure it runs within a reasonable timeframe.
Maximization: Increase the number of usable primers stored within a DNA pool.
Algorithm Profiling: Analyze and profile the algorithm to understand its bottlenecks and scalability.
Approach


**Naive Algorithm**: Started with a basic Python implementation that was tested to generate primers but proved too slow for practical use.
Bioinformatics Tools:

##### **Primer3**: Used for melting temperature calculation.

##### **NUPACK**: For secondary structure validation.

##### **CD-HIT**: For sequence similarity filtering.

### Advanced Strategy

Primer Generation: N random primers are generated, followed by filtering against specific intra-primer restrictions
Inter-Primer Checks: The algorithm uses a trie-based system to efficiently filter primers based on Hamming distance, avoiding unnecessary comparisons.
and for the inter-complementarity , by utilizing sets to store reverse complement patterns, the solution achieves faster inter-primer validation.

Key Biological Constraints

###### GC Content: Between 45% and 55%.

###### Long Homo-polymers: Disallowed to prevent reading errors.

###### Sequence Complementarity: Both intra-primer and inter-primer self-complementarity are minimized to avoid unintended interactions.

###### Hamming Distance: Ensures primers are sufficiently different to avoid misidentification.

Tools & Technologies

This approach effectively scales for large datasets,The algorithm performs significantly faster compared to basic methods, especially in filtering primers for reverse complementarity and Hamming distance checks.

## Contributors

##### Saar Cohen

##### Shahar Trabelsi

#### Advisor: Dr. Sarel Cohen

#### Collaboration: Prof. Dalit Naor

References

Lee Organick, et al., "Random Access in Large-Scale DNA Data Storage", Nature Biotechnology, 2018.
Additional bioinformatics resources: DNA Storage Alliance


results:
![final_results_1and10mil.jpg](extra%2Ffinal_results_1and10mil.jpg)
![Filtration Pipeline.jpg](extra%2FFiltration%20Pipeline.jpg)
![sorted_unsorted.jpg](extra%2Fsorted_unsorted.jpg)
![timing.jpg](extra%2Ftiming.jpg)
![hamming_distance with trie-tree.jpg](extra%2Fhamming_distance%20with%20trie-tree.jpg)