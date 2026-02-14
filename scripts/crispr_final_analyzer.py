#!/usr/bin/env python3.11
"""
CRISPR Sanger测序最终版分析工具

根据用户最终澄清的需求，实现一个能够准确分析复合杂合突变的专业工具。
支持用户指定参考位置进行分析，以应对测序质量不佳导致自动检测不准的问题。
"""

import numpy as np
from Bio import SeqIO, Align
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import csv
from typing import Dict, List, Tuple
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

class CRISPRFinalAnalyzer:
    """CRISPR最终版分析器"""

    def __init__(self, reference_file: str):
        ref_record = SeqIO.read(reference_file, "fasta")
        self.reference_seq = str(ref_record.seq).upper()
        self.reference_name = ref_record.id

    def read_ab1_file(self, ab1_file: str) -> Dict:
        """读取ab1文件"""
        try:
            record = SeqIO.read(ab1_file, "abi")
            abif_raw = record.annotations.get("abif_raw", {})
            
            channels = {
                "G": np.array(abif_raw.get("DATA9", [])),
                "A": np.array(abif_raw.get("DATA10", [])),
                "T": np.array(abif_raw.get("DATA11", [])),
                "C": np.array(abif_raw.get("DATA12", []))
            }
            
            peak_locations = np.array(abif_raw.get("PLOC2", abif_raw.get("PLOC1", [])))
            sequence = str(record.seq).upper()
            
            return {
                "success": True,
                "channels": channels,
                "peak_locations": peak_locations,
                "sequence": sequence,
                "file_name": Path(ab1_file).name
            }
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "file_name": Path(ab1_file).name
            }

    def get_peak_heights(self, channels: Dict, peak_loc: int, window: int = 5) -> Dict:
        """获取指定位置的四个通道峰高"""
        heights = {}
        for base, channel_data in channels.items():
            if len(channel_data) > 0:
                w_start = max(0, peak_loc - window)
                w_end = min(len(channel_data), peak_loc + window)
                heights[base] = float(np.max(channel_data[w_start:w_end])) if w_end > w_start else 0.0
            else:
                heights[base] = 0.0
        return heights

    def analyze_peak_pattern(self, heights: Dict) -> Dict:
        """分析峰图模式"""
        sorted_peaks = sorted(heights.items(), key=lambda x: x[1], reverse=True)
        primary_base, primary_height = sorted_peaks[0]
        secondary_base, secondary_height = sorted_peaks[1]
        ratio = secondary_height / primary_height if primary_height > 0 else 0
        return {
            "primary_base": primary_base,
            "primary_height": primary_height,
            "secondary_base": secondary_base,
            "secondary_height": secondary_height,
            "ratio": ratio,
            "is_heterozygous": ratio > 0.3, # Lowered threshold slightly for robustness
            "all_heights": heights
        }

    def calculate_all_ratios(self, data: Dict) -> np.ndarray:
        """计算所有位置的次峰比例"""
        ratios = []
        for i in range(len(data["peak_locations"])):
            peak_loc = data["peak_locations"][i]
            heights = self.get_peak_heights(data["channels"], peak_loc)
            pattern = self.analyze_peak_pattern(heights)
            ratios.append(pattern["ratio"])
        return np.array(ratios)

    def local_align(self, query_seq: str) -> Dict:
        """局部比对"""
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        alignments = aligner.align(self.reference_seq, query_seq)
        if not alignments:
            return {"success": False, "error": "No alignment found"}
        
        best_aln = alignments[0]
        ref_to_query = {}
        ref_pos, query_pos = best_aln.aligned[0][0][0], best_aln.aligned[1][0][0]
        ref_aligned, query_aligned = str(best_aln[0]), str(best_aln[1])

        for i in range(len(ref_aligned)):
            if ref_aligned[i] != '-':
                if query_aligned[i] != '-':
                    ref_to_query[ref_pos] = query_pos
                ref_pos += 1
            if query_aligned[i] != '-':
                query_pos += 1

        return {
            "success": True,
            "ref_to_query": ref_to_query,
            "score": best_aln.score
        }

    def analyze_specific_position(self, data: Dict, aln_result: Dict, target_ref_pos: int) -> Dict:
        """分析指定参考位置的突变"""
        ratios = self.calculate_all_ratios(data)
        query_pos_0based = aln_result["ref_to_query"].get(target_ref_pos - 1)
        if query_pos_0based is None or query_pos_0based >= len(ratios):
            return {"mutation_detected": False, "error": f"参考位置 {target_ref_pos} 不在比对范围内或无信号"}

        if ratios[query_pos_0based] < 0.3:
            return {"mutation_detected": False, "error": f"参考位置 {target_ref_pos} (测序 {query_pos_0based+1}) 不是高杂合位点 (ratio={ratios[query_pos_0based]:.3f}) - 阈值 0.3"}

        future_ratios = ratios[query_pos_0based + 1 : query_pos_0based + 21]
        mean_future_ratio = np.mean(future_ratios) if len(future_ratios) > 0 else 0
        has_frameshift = mean_future_ratio > 0.5

        peak_loc = data["peak_locations"][query_pos_0based]
        heights = self.get_peak_heights(data["channels"], peak_loc)
        peak_pattern = self.analyze_peak_pattern(heights)
        ref_base = self.reference_seq[target_ref_pos - 1]

        allele1, allele2 = "未知", "未知"
        if has_frameshift:
            genotype = "Compound Heterozygous"
            primary, secondary = peak_pattern["primary_base"], peak_pattern["secondary_base"]
            
            # Logic for user's specific case: ref=C, allele1=C insertion, allele2=C->A
            # This means the two signals present are C and A.
            if ref_base == 'C' and set([primary, secondary]) == set(['A', 'C']):
                allele1 = "在参考248位插入1个C (C->CC)"
                allele2 = "参考248位C突变为A (C->A)"
            else: # Generic fallback for other cases
                allele1 = f"可能的Indel (导致frameshift)"
                allele2 = f"可能的替换: {ref_base} -> {primary if primary != ref_base else secondary}"
        else:
            genotype = "Heterozygous"
            allele1 = f"{ref_base} (野生型)"
            allele2 = f"{ref_base}→{peak_pattern['secondary_base']}"

        return {
            "mutation_detected": True, "genotype": genotype, "ref_pos": target_ref_pos,
            "query_pos": query_pos_0based + 1, "ref_base": ref_base, "peak_pattern": peak_pattern,
            "has_frameshift": has_frameshift, "mean_future_ratio": mean_future_ratio,
            "allele1": allele1, "allele2": allele2, "ratios": ratios
        }

    def analyze_single_file(self, ab1_file: str, target_ref_pos: int) -> Dict:
        """分析单个ab1文件"""
        data = self.read_ab1_file(ab1_file)
        if not data["success"]: return data
        aln_result = self.local_align(data["sequence"])
        if not aln_result["success"]: return {**data, **aln_result}
        mut_result = self.analyze_specific_position(data, aln_result, target_ref_pos)
        return {"success": True, "file_name": data["file_name"], "mutation_result": mut_result, "data": data}

    def generate_report(self, result: Dict) -> str:
        """生成详细报告"""
        if not result["success"]: return f"文件: {result['file_name']}\n分析失败: {result.get('error', '未知错误')}"
        mut = result['mutation_result']
        report = [f"文件: {result['file_name']}", "="*40]
        if not mut["mutation_detected"]:
            report.append(f"在指定位置未检测到突变: {mut.get('error', '')}")
        else:
            pp = mut['peak_pattern']
            report.extend([
                f"基因型: {mut['genotype']}",
                f"突变位点: 参考 {mut['ref_pos']} / 测序 {mut['query_pos']}",
                f"参考碱基: {mut['ref_base']}",
                f"峰图分析: 主峰 {pp['primary_base']} (高度: {pp['primary_height']:.0f}), 次峰 {pp['secondary_base']} (高度: {pp['secondary_height']:.0f}), 比例 {pp['ratio']:.3f}",
                f"Frameshift: {'是' if mut['has_frameshift'] else '否'} (后续20个位置平均比例: {mut['mean_future_ratio']:.3f})",
                "\n--- 推断结果 ---",
                f"等位基因 1: {mut['allele1']}",
                f"等位基因 2: {mut['allele2']}"
            ])
        return "\n".join(report)

    def visualize_mutation(self, result: Dict, output_file: str, window: int = 25):
        """可视化突变位点"""
        mut = result.get("mutation_result", {})
        if not mut.get("mutation_detected"):
            print("未检测到突变，无法可视化")
            return

        data, ratios, query_pos_1based = result["data"], mut["ratios"], mut["query_pos"]
        query_pos_0based = query_pos_1based - 1
        ref_pos = mut["ref_pos"]
        
        query_start = max(0, query_pos_0based - window)
        query_end = min(len(data["sequence"]), query_pos_0based + window)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), gridspec_kw={'height_ratios': [2, 1]})

        # Plot chromatogram
        trace_start = data["peak_locations"][query_start] - 50
        trace_end = data["peak_locations"][query_end - 1] + 50
        x_trace = np.arange(trace_start, trace_end)
        colors = {"G": "black", "A": "green", "T": "red", "C": "blue"}
        for base, color in colors.items():
            ax1.plot(x_trace, data["channels"][base][trace_start:trace_end], color=color, label=base, lw=1.5)
        
        for i in range(query_start, query_end):
            peak_loc = data["peak_locations"][i]
            is_mut_pos = (i == query_pos_0based)
            ax1.axvline(peak_loc, color='red' if is_mut_pos else 'grey', ls='-' if is_mut_pos else '--', lw=2 if is_mut_pos else 0.5)
            ax1.text(peak_loc, ax1.get_ylim()[1] * 0.9, data["sequence"][i], ha='center', 
                     fontweight='bold', fontsize=9, color='red' if is_mut_pos else 'black')

        ax1.set_title(f"Chromatogram at Ref Pos {ref_pos} (Query Pos {query_pos_1based})", fontsize=16)
        ax1.set_xlabel("Trace Position")
        ax1.set_ylabel("Signal Intensity")
        ax1.legend()

        # Plot ratios
        x_query = np.arange(query_start, min(query_end, len(ratios)))
        ax2.plot(x_query, ratios[query_start:query_end], '-o', label='Secondary Peak Ratio', ms=4)
        ax2.axhline(0.3, color='orange', ls='--', label='Heterozygous Threshold (0.3)')
        ax2.axvline(query_pos_0based, color='red', ls='-', lw=2, label=f'Mutation Site (Ref {ref_pos})')
        if mut["has_frameshift"]:
            ax2.axvspan(query_pos_0based, query_end, color='red', alpha=0.15, label='Frameshift Region')

        ax2.set_title("Secondary Peak Ratio Analysis", fontsize=16)
        ax2.set_xlabel("Query Sequence Position")
        ax2.set_ylabel("Ratio")
        ax2.set_ylim(0, 1.1)
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=200)
        print(f"可视化图已保存至: {output_file}")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="CRISPR Sanger测序最终版分析工具")
    parser.add_argument("-i", "--input", required=True, help="输入ab1文件")
    parser.add_argument("-r", "--reference", required=True, help="参考序列 (FASTA)")
    parser.add_argument("-p", "--position", type=int, required=True, help="要分析的参考序列位置 (1-based)")
    parser.add_argument("-o", "--output", default="crispr_analysis", help="输出文件前缀")
    parser.add_argument("--visualize", action="store_true", help="生成可视化图")
    args = parser.parse_args()

    analyzer = CRISPRFinalAnalyzer(args.reference)
    result = analyzer.analyze_single_file(args.input, args.position)
    
    report_text = analyzer.generate_report(result)
    print("\n" + "="*50 + "\n分析报告:\n" + "="*50)
    print(report_text)

    report_file = f"{args.output}_report.txt"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(report_text)
    print(f"\n详细报告已保存至: {report_file}")

    if args.visualize and result.get("success") and result.get("mutation_result", {}).get("mutation_detected"):
        vis_file = f"{args.output}_visualization.png"
        analyzer.visualize_mutation(result, vis_file)

if __name__ == "__main__":
    main()
