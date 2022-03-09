<h2 class="c4" id="h.70zate8h3f8a"><span>Title: </span><span>RNA-</span><span>Health</span><sup><a href="#cmnt1"
            id="cmnt_ref1">[a]</a></sup></h2>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3"><span class="c14">Github handle</span><span>: </span><span class="c9"><a class="c10"
            href="https://www.google.com/url?q=http://github.com/antonkulaga/rna-seq&amp;sa=D&amp;source=editors&amp;ust=1646824004878760&amp;usg=AOvVaw1DG1c6VgsabEKj-4kgDZ9D">http://github.com/antonkulaga/rna-seq</a></span>
</p>
<h3 class="c12" id="h.ous3fw7b4a0u"><span class="c19">Who are we?</span></h3>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c20"><span>Our team is formed of a mix of multi-disciplinary scientists with strong expertise in
        bioinformatics, systems biology, genetics and molecular cell biology. The group is involved in both academic and
        industry activities. Academic projects in the field of aging-research are generally carried out by our team as
        part of the Systems Biology of Aging Group (</span><span class="c9"><a class="c10"
            href="https://www.google.com/url?q=http://aging-research.group&amp;sa=D&amp;source=editors&amp;ust=1646824004879701&amp;usg=AOvVaw1T1wqT79j1LSldXhdt3zXG">www.aging-research.group</a></span><span>)
        at the Institute of Biochemistry of the Romanian Academy. More entrepreneurial projects are implemented within a
        spin-off company, CellFabrik (</span><span class="c9"><a class="c10"
            href="https://www.google.com/url?q=http://www.cellfabrik.bio&amp;sa=D&amp;source=editors&amp;ust=1646824004879888&amp;usg=AOvVaw0q9_g_qoNIL2PTvWu_SJli">www.cellfabrik.bio</a></span><span
        class="c2">). &nbsp;</span></p>
<p class="c3"><span class="c2">The group&#39;s expertise includes analysis of transcriptomics [1][2], single-cell RNA
        sequencing [2], machine learning [1][3], and aging biology[4][5].</span></p>
<p class="c3 c8"><span class="c2"></span></p>
<h3 class="c12" id="h.y463f3t6rjs8"><span class="c19">What is RNA-Health</span></h3>
<h4 class="c24" id="h.gh96o5r216zy"><span class="c26">Sequencing and health</span></h4>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3"><span class="c2">In biology there are multiple types of sequencing that reveal different health aspects,
        for example:</span></p>
<ul class="c1 lst-kix_a8g00ndx3zsl-0 start">
    <li class="c0 li-bullet-0"><span class="c2">Genomics sequencing can tell about your hereditary health and some
            risks. </span></li>
    <li class="c0 li-bullet-0"><span class="c2">Methylation (bisulfite) sequencing shows epigenetic states of the cells
            and can also be useful for measuring your biological age</span></li>
    <li class="c0 li-bullet-0"><span class="c2">Adaptive immune repertoire sequencing (AIRR-Seq) can show which
            antibodies and t-cell receptors your immune cells produce in response to bacterial, viral, autologous and
            other antigens</span></li>
    <li class="c0 li-bullet-0"><span class="c2">Proteomic assays (SomaScan and others) show current state of some organs
            but provide no hereditary info</span></li>
    <li class="c0 li-bullet-0"><span class="c2">Etc.</span></li>
</ul>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3 c8"><span class="c2"></span></p>
<h3 class="c12" id="h.rmac2rb6mlbw"><span class="c19">RNA sequencing</span></h3>
<p class="c3"><span
        style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 718.82px; height: 608.00px;"><img
            alt="" src="images/image1.png"
            style="width: 718.82px; height: 608.00px; margin-left: 0.00px; margin-top: 0.00px; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px);"
            title=""></span></p>
<p class="c20"><span class="c2">In comparison with other types of sequencing, RNA-Seq is the &quot;Jack of all trades
        but master of none&quot;. It can partially substitute other sequencing types:</span></p>
<ul class="c1 lst-kix_z5ii4dgf0cw-0 start">
    <li class="c5 li-bullet-0"><span>RNA variant-calling can tell about coding genetic variations in your genome
            [6]</span><sup><a href="#cmnt4" id="cmnt_ref4">[d]</a></sup></li>
    <li class="c5 li-bullet-0"><span>RNA quantification can tell about the protein content and cell types proportions
            (cell types deconvolution) in the blood and can also be used to compare different conditions (differential
            expression)</span><sup><a href="#cmnt5" id="cmnt_ref5">[e]</a></sup></li>
    <li class="c5 li-bullet-0"><span>Adaptive immune repertoire (antibodies and T cell receptors) can be reconstructed
            from a subset of RNA [7]</span><sup><a href="#cmnt6" id="cmnt_ref6">[f]</a></sup></li>
    <li class="c5 li-bullet-0"><span>Aging clock can also be computed from gene expressions [8]</span><sup><a
                href="#cmnt7" id="cmnt_ref7">[g]</a></sup></li>
</ul>
<p class="c20"><span class="c2">In the same time, doing multiple sequencing on yourself or your pets can be
        prohibitively expensive while blood RNA-Seq is relatively cheap and can be useful for multiple purposes and can
        also integrate with other types of sequencing. This makes RNA-Seq promising tool for quick health exploration or
        refinement of information from other sequencing modalities.</span></p>
<p class="c20 c8"><span class="c2"></span></p>
<h3 class="c18" id="h.23yr2g3o66la"><span class="c19">RNA-Health</span></h3>
<p class="c20 c8"><span class="c2"></span></p>
<p class="c20"><span class="c2">The goal of this project is to provide an opensource toolbox that will allow:</span></p>
<ul class="c1 lst-kix_f5m5nbtoy9fy-0 start">
    <li class="c5 li-bullet-0"><span>RNA age prediction [8] which </span><span>involves training RNA</span><sup><a
                href="#cmnt8" id="cmnt_ref8">[h]</a></sup><sup><a href="#cmnt9" id="cmnt_ref9">[i]</a></sup><sup><a
                href="#cmnt10" id="cmnt_ref10">[j]</a></sup><span class="c2">&nbsp;aging clock on available public
            datasets as well as integrating known clocks like BitAge [8]</span></li>
    <li class="c5 li-bullet-0"><span class="c2">Infering coding genes variations relevant for the health and providing
            polygenic scores for diseases (collaboration with Just-DNA-Seq project)</span></li>
    <li class="c5 li-bullet-0"><span class="c2">Extraction of adaptive immune repertoires from RNA-Seq data and
            generating adaptive immune report. Processing of public RNA-Seq data of different</span></li>
    <li class="c5 li-bullet-0"><span class="c2">Refinement of genetic variations extracted from DNA-Seq with RNA-Seq
            data (collaboration with Just-DNA-Seq project)</span></li>
    <li class="c5 li-bullet-0"><span class="c2">Privacy-preserving storage and analysis of your personal RNA-Seq
            (collaboration with Genomes-DAO)</span></li>
    <li class="c5 li-bullet-0"><span>Integrating the pipelines with Lab DAO ecosystem that will allow everybody run them
            online</span></li>
</ul>
<p class="c20 c8"><span class="c2"></span></p>
<h3 class="c12" id="h.qks8il8mpyi1"><span class="c19">Why a Gitcoin Grant?</span></h3>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c20"><span class="c2">This Gitcoin grant will be used to fund and speed up the development of the project.
        Often, when research code is written there is a lack of incentive to turn it from just an executable supplement
        of the academic papers to a library and a toolbox that everybody can use benefit from. Gitcoin funding will
        allow us to turn our previous research pipelines in an extensible opensource toolbox that other researchers,
        citizen scientist and quantify-self enthusiasts can use, extend and apply both to their personal data or
        datasets of their interest.</span></p>
<p class="c3 c8"><span class="c11"></span></p>
<h3 class="c12" id="h.p4fcf2hii8bw"><span>Other ways to help us</span></h3>
<ul class="c1 lst-kix_lhh5fod2035p-0 start">
    <li class="c0 li-bullet-0"><span class="c2">If you have programming skills you can help with the development</span>
    </li>
    <li class="c0 li-bullet-0"><span class="c2">If you have biological skills you can help with prioritising genes in
            genetic reports and literature curation</span></li>
    <li class="c0 li-bullet-0"><span>If you are analysing your data you can always provide feedback on what should be
            improved</span></li>
</ul>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3 c8"><span class="c2"></span></p>
<p class="c3"><span class="c2">References:</span></p>
<p class="c3 c8"><span class="c2"></span></p>
<ol class="c1 lst-kix_5e3ecna6101l-0 start" start="1">
    <li class="c0 li-bullet-0"><span class="c2">AY Kulaga, E Ursu, D Toren, V Tyshchenko, R Guinea&hellip; -
            International journal of molecular sciences, 2021</span></li>
    <li class="c0 li-bullet-0"><span>Lagger, C., Ursu, E., Equey, A., Avelar, R. A., Pisco, A. O., Tacutu, R., &amp; de
            Magalh&atilde;es, J. P. (2021). scAgeCom: a murine atlas of age-related changes in intercellular
            communication inferred with the package scDiffCom. </span><span class="c6">bioRxiv</span><span
            class="c2">.</span></li>
    <li class="c0 li-bullet-0"><span class="c9"><a class="c10"
                href="https://www.google.com/url?q=https://scholar.google.ro/scholar?oi%3Dbibs%26cluster%3D6402377544428486051%26btnI%3D1%26hl%3Den&amp;sa=D&amp;source=editors&amp;ust=1646824004884765&amp;usg=AOvVaw3MMo1rzpeeX0ZsGoPb9v9Z">Learning
                flat representations with artificial neural networks</a></span><span>&nbsp;</span><span class="c2">V
            Constantinescu, C Chiru, T Boloni, A Florea&hellip; - Applied Intelligence, 2021</span></li>
    <li class="c0 li-bullet-0"><span>Bunu, Gabriela, et al. &quot;SynergyAge, a curated database for synergistic and
            antagonistic interactions of longevity-associated genes.&quot; </span><span class="c6">Scientific
            data</span><span class="c2">&nbsp;7.1 (2020): 1-11.</span></li>
    <li class="c0 li-bullet-0"><span>Tacutu, Robi, et al. &quot;Human ageing genomic resources: new and updated
            databases.&quot; </span><span class="c6">Nucleic acids research</span><span class="c2">&nbsp;46.D1 (2018):
            D1083-D1090.</span></li>
    <li class="c0 li-bullet-0"><span>Koboldt, Daniel C. &quot;Best practices for variant calling in clinical
            sequencing.&quot; </span><span class="c6">Genome Medicine</span><span class="c2">&nbsp;12.1 (2020):
            1-13.</span></li>
    <li class="c0 li-bullet-0"><span>Song, Li, et al. &quot;TRUST4: immune repertoire reconstruction from bulk and
            single-cell RNA-seq data.&quot; </span><span class="c6">Nature Methods</span><span>&nbsp;18.6 (2021):
            627-630.</span></li>
    <li class="c0 li-bullet-0"><span>Meyer, David H., and Bj&ouml;rn Schumacher. &quot;BiT age: A
            transcriptome&#8208;based aging clock near the theoretical limit of accuracy.&quot; </span><span
            class="c6">Aging cell</span><span>&nbsp;20.3 (2021): e13320.</span></li>
</ol>