$(function() {
  for (cat in data) {
    var html = "";
    for (i = 0; i < data[cat].length; i++)
      if (cat != "cancers")
        html += `<option value="${data[cat][i]}">${data[cat][i]}</option>`;
      else
        html += `<option value="${data[cat][i]}">TCGA-${data[cat][i]}</option>`;

    
    $(`.${cat}-select`).append(html);
    $(`.${cat}-select`).selectpicker("render");
    $(`.${cat}-spinner`).remove();
  }
  
  $(".cancer-search .btn").click(function() {
    var cancer = $(".cancer-search select").val();
    if (cancer == "")
      alert("Please select a cancer");
    else
      window.location.href = `/cancer/${cancer}`;
  });

  $(".miRNA-search .btn").click(function() {
    var miRNA = $(".miRNA-search select").val();
    if (miRNA == "")
      alert("Please select a miRNA");
    else
      window.location.href = `/miRNA/${miRNA}`;
  });
  
  $(".isomiR-search .btn").click(function() {
    var isomiR = $(".isomiR-search select.isomiRs-select").val();
    var cancer = $(".isomiR-search select.cancers-select").val();
    if ((isomiR == "") || (cancer == "")) {
      alert("Please select an isomiR and a cancer");
    } else {
      if (cancer == "pan-cancer")
        window.location.href = `/molecule/${isomiR}`;
      else
        window.location.href = `/cancer_molecule/${cancer}/${isomiR}`;
    }
  });
  
  $(".gene-search .btn").click(function() {
    var gene = $(".gene-search select.genes-select").val();
    var cancer = $(".gene-search select.cancers-select").val();
    if ((gene == "") || (cancer == "")) {
      alert("Please select a gene and a cancer");
    } else {
      if (cancer == "pan-cancer")
        window.location.href = `/molecule/${gene}`;
      else
        window.location.href = `/cancer_molecule/${cancer}/${gene}`;
    }
  });
  
  $(".triple-search .btn").click(function() {
    var isomiR = $(".triple-search select.isomiRs-select").val();
    var gene = $(".triple-search select.genes-select").val();
    var cancer = $(".triple-search select.cancers-select").val();
    if ((isomiR == "") || (gene == "") || (cancer == "")) 
      alert("Please select an isomiR, a gene and a cancer");
    else
      window.location.href = `/cancer_isomir_target/${cancer}/${isomiR}/${gene}`;
  });
});
